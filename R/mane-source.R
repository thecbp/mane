IND = function(n_subj, n_trts, n_periods, n_obvs,
               betas, y_sigma, chains, iter, lag,
               stabilize = NULL,
               objective = "max",
               adapt_delta = 0.999,
               max_treedepth = 15,
               seed = 1) {

  set.seed(seed)

  # Stopping the RMD CHECK for global variables
  id = NULL

  # Storing the posterior samples for after each update
  trial_posteriors = list()

  print("Burn in period")

  # Initial FRN phase data for burn-in
  start_data = generate_FRN_data(n_subj, n_trts, n_obvs, betas, y_sigma)

  # Build the formula for the model (X1 used as reference)
  main = paste0("X", 2:n_trts, collapse = " + ")
  f = paste0("Y ~ ", main, " + ar(time = time, p = ", lag, ")")

  # Setting priors for the model
  individual_priors = c(
    brms::set_prior("normal(0, 100)", class = "Intercept"),
    brms::set_prior("normal(0, 100)", class = "b"),
    brms::set_prior("normal(0, 100)", class = "ar")
    # Using the default priors for sigma (half_t)
  )

  # Choosing to use the default priors for the correlation of random effects

  # Generating posterior samples after data collected for everyone in period
  posteriors = start_data %>%
    tidyr::nest(data = colnames(start_data)[-1]) %>% # group data on id column
    dplyr::mutate(
      model = map(data, function(df) {
        brms::brm(data = df, formula = f,
                  family = gaussian(link = "identity"),
                  chains = chains, iter = iter,
                  prior = individual_priors,
                  control = list(max_treedepth = max_treedepth,
                                 adapt_delta = adapt_delta),
                  refresh = 0)

      })
    )

  # Record all of the models
  trial_posteriors[[n_trts]] = posteriors

  # If stabilize not initalized, use Thall & Walthen's optimal parameter
  if (is.null(stabilize)) { c = n_trts / (2 * n_periods) }
  else { c = stabilize }

  # Calculate allocation probabilities based on Thompson Sampling
  current_probs = allocate_probabilities(posteriors, c, objective = objective)
  current_probs$period = n_trts # By convention, just assign the last burn in period

  # Start collecting all of the allocation probabilities into one place
  trial_probs = current_probs

  # Index to start the loop, since n_trts periods have been used in FRN phase
  adaptive_start = n_trts + 1 # FRN requires one period from every treatment

  # For each treatment cycle...
  for (per in adaptive_start:n_periods) {

    print(paste0("Starting period ", per))

    # Randomize each subject to new treatment, sample data, add to current set
    updated_posteriors = posteriors %>%
      dplyr::mutate(
        next_trt = purrr::map_dbl(id, function(i) {

          id_probs = current_probs %>%
            dplyr::filter(id == i) %>%
            dplyr::select(tidyselect::starts_with("X")) %>%
            unlist()

          sample(1:n_trts, size = 1, prob = id_probs)

        }),
        newdata = purrr::map2(next_trt, id, function(nt, i) {

          x = matrix(0, nrow = 1, ncol = n_trts)
          x[1] = 1        # intercept
          x[nt] = 1       # setting next treatment
          design = x[rep(1, n_obvs), ]
          id_obvs = rep(x %*% betas[i, ], n_obvs)
          id_y = id_obvs + stats::rnorm(n_obvs, 0, y_sigma)

          if (n_obvs == 1) {
            id_data = tibble::as_tibble(
              list(
                Y = id_y,
                period = per
              ))
            for (j in seq_len(n_trts)) {
              id_data[[paste0("X", j)]] = x[j]
            }
          } else {
            id_data = cbind(design, id_y, rep(per, n_obvs)) %>% tibble::as_tibble()
            colnames(id_data) = c(paste0("X", 1:n_trts), "Y", "period")
          }

          id_data

        }),
        full = purrr::map2(data, newdata, function(d, nd) {

          # Error has something to do with data
          full_data = dplyr::bind_rows(d, nd)
          full_data$time = 1:nrow(full_data)
          full_data

        }),
        model = purrr::map(full, function(df) {

          brms::brm(data = df, formula = f,
                    family = gaussian(link = "identity"),
                    chains = chains, iter = iter,
                    prior = individual_priors,
                    control = list(max_treedepth = max_treedepth,
                                   adapt_delta = adapt_delta),
                    refresh = 0)

        })
      ) %>%
      dplyr::select(id, data = full, model)

    # Record and update the posteriors used in the loop
    trial_posteriors[[per]] = updated_posteriors
    posteriors = updated_posteriors

    # If stabilize not initialized, use Thall & Walthen's optimal parameter
    if (is.null(stabilize)) { c = per / (2 * n_periods) }
    else { c = stabilize }

    # Readjust reallocation probabilities based on new updated posterior
    current_probs = allocate_probabilities(posteriors, c, objective = objective)
    current_probs$period = per

    # Record the new probabilities
    trial_probs = dplyr::bind_rows(trial_probs, current_probs)

  } # end of period loop

  # Set up output object
  output = list(
    trial_params = list(
      n_subj = n_subj,
      n_periods = n_periods,
      n_obvs = n_obvs,
      betas = betas,
      y_sigma = y_sigma
    ),
    stan_params = list(
      chains = chains,
      iter = iter,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    ),
    allocation = trial_probs,
    posteriors = trial_posteriors
  )

  output

}

generate_FRN_data = function(n_subj,
                             n_trts,
                             n_obvs,
                             betas,
                             y_sigma) {

  # Matrix of single person getting all of the treatments for n_obvs each
  trt = diag(n_trts)
  trt[, 1] = 1
  colnames(trt) = paste0("X", 1:n_trts)
  trt = trt[rep(seq_len(nrow(trt)), each = n_obvs),] # single design matrix
  X = trt[rep(seq_len(nrow(trt)), each = n_subj),] # matrices for each person

  # Calculate outcome with given betas (column is person, row is outcome for first obvs)
  # Also adding within-subject noise
  Y = ((trt %*% t(betas)) %>% c()) + stats::rnorm(nrow(X), 0, y_sigma)

  id = rep(1:n_subj, times = n_obvs * n_trts)

  output = cbind(id = id, X, Y = Y) %>% tibble::as_tibble()
  output = output %>%
    dplyr::mutate(
      period = rep(1:n_trts, each = n_subj * n_obvs)
    ) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(
      time = 1:(n_trts * n_obvs) # Mark the time of observation
    ) %>%
    ungroup()

  output

}

allocate_probabilities = function(posterior, c, objective = "max") {

  # Set the correct function to maximize on
  if (objective == "max") { m = max }
  else { m = min }

  output = posterior %>%
    dplyr::transmute(
      id = id,
      probs = purrr::map(model, function(mod) {

        betas = as.data.frame(mod) %>% dplyr::select(tidyselect::contains("b_"))
        n_trts = ncol(betas)

        # Design matrix representing outcome based on treatments
        trt = diag(n_trts)
        trt[, 1] = 1
        colnames(trt) = paste0("X", 1:n_trts)

        # Find which arm produced best treatment effects
        posterior_outcomes = trt %*% t(betas)
        find_optimal = function(x) {which(x == m(x))}
        optimal_arms = apply(posterior_outcomes, 2, find_optimal)

        probs = c(table(optimal_arms) / nrow(betas))

        # Stabilize the probabilities
        stab_probs = probs^c / sum(probs^c)

        prob_df = stab_probs %>%
          tibble::as_tibble() %>%
          dplyr::mutate(
            trt = names(stab_probs)
          ) %>%
          tidyr::pivot_wider(names_from = trt, values_from = value)

        prob_df
      })
    ) %>%
    tidyr::unnest(probs)

  # Adjust column names
  n_trts = ncol(output) - 1
  colnames(output) = c("id", paste0("X", 1:n_trts))

  # Replace NA with 0
  replace_list = list()
  for (col in paste0("X", 1:n_trts)) {
    replace_list[[col]] = 0
  }

  output = output %>% replace_na(replace_list)

  output
}
