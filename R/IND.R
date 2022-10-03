#' Simulate an individual bandit design for a multi-arm N-of-1 trial
#'
#' `IND()` simulates an individual bandit design for a multi-arm N-of-1 trial.
#' This design comes in two phases: a burn-in phase where all the treatments are
#' gievn in a random order, and an adaptive phase where the next treatment
#' randomized is chosen via Thompson Sampling. The allocation probabilities are
#' calculated via the posterior probability of each treatment being optimal.
#'
#' In the individual bandit design, the posterior treatment effects are
#' updated after all individual have finished a treatment period. It is assumed
#' that estimating the individual-level effects are of interest.
#'
#' @inheritParams FRN
#'
#' @return List containing the simulation parameters and resulting trial data that came from the parameters
#' @export
#'
#' @examples
#' \dontrun{
#' # Generating data for simulation
#' set.seed(1)
#'
#' # Covariance between betas, factorized into scalar and correlation components
#' tau = 2 * diag(3)
#' Omega = diag(3)
#' Omega[1,2] = Omega[2,1] = 0.1
#' Omega[1,3] = Omega[3,1] = 0.1
#' Omega[3,2] = Omega[2,3] = 0.1
#' Sigma = tau %*% Omega %*% tau
#'
#' # Generating individual treatment effects from unique distributions
#' betas = array(0, dim = c(2, 3))
#' betas[1,] = MASS::mvrnorm(n = 1, mu = c(3, -2, 4), Sigma = Sigma)
#' betas[2,] = MASS::mvrnorm(n = 1, mu = c(1, 0, 3), Sigma = Sigma)
#' betas = round(betas, 1)
#' stanfile = 'stan/model_ind.stan'
#'
#' # ind = IND(n_subj = 2, n_periods = 6, n_obvs = 5, n_trts = 3, y_sigma = 2,
#' #           stanfile = stanfile, betas = betas, chains = 1, warmup = 1000,
#' #           iter = 3000)
#' }
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
    tidyr::nest(data = colnames(current_data)[-1]) %>% # group data on id column
    dplyr::mutate(
      model = map(data, function(df) {
        brms::brm(data = df, formula = f,
                  family = gaussian(link = "identity"),
                  chains = chains, iter = iter,
                  prior = individual_priors,
                  control = list(max_treedepth = max_treedepth,
                                 adapt_delta = adapt_delta))
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
                                   adapt_delta = adapt_delta))

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


# Update the posterior distribution after observing these new values
# posterior = rstanarm::stan_glmer(
#   formula = model_formula,
#   data = current_data, family = "gaussian",
#   prior_intercept = rstanarm::normal(100, 10),
#   prior = rstanarm::normal(0, 2.5),
#   prior_aux = rstanarm::exponential(1, autoscale = TRUE),
#   prior_covariance = rstanarm::decov(reg = 1, conc = 1, shape = 1, scale = 1),
#   chains = chains, iter = iter, seed = seed, adapt_delta = adapt_delta,
#   control = list(max_treedepth = max_treedepth)
# )
Footer
