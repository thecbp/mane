#' Simulate an aggregated bandit design for a multi-arm N-of-1 trial
#'
#' `AGG()` simulates an aggregated bandit design for a multi-arm N-of-1 trial.
#' This design comes in two phases: a burn-in phase where all the treatments are
#' gievn in a random order, and an adaptive phase where the next treatment
#' randomized is chosen via Thompson Sampling. The allocation probabilities are
#' calculated via the posterior probability of each treatment being optimal.
#'
#' In an aggregated bandit design, the posterior treatment effects are updated
#' after an individual finishes a treatment period. It is assumed that
#' estimating the population-level effects are of interest.
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
#' # Generating individual treatment effects from population level
#' betas = MASS::mvrnorm(n = 2, mu = c(3, -2, 4), Sigma = Sigma)
#' betas = round(betas, 1)
#' stanfile = 'stan/model_agg.stan'
#'
#' # agg = AGG(n_subj = 2, n_periods = 6, n_obvs = 5, n_trts = 3, y_sigma = 2,
#' #           stanfile = stanfile, betas = betas, chains = 1, warmup = 1000,
#' #           iter = 3000)
#' }
AGG = function(n_subj, n_trts, n_periods, n_obvs, betas, y_sigma, optimize = "max",
               chains,  iter, adapt_delta = 0.999, max_treedepth = 17, seed = 1) {

  # Stopping the RMD CHECK for global variables
  id = NULL

  # Storing the posterior samples for after each update
  trial_posteriors = list()

  print("Burn in period")
  # Initial FRN phase data for burn-in
  current_data = generate_FRN_data(n_subj, n_trts, n_periods = 1, n_obvs, betas, y_sigma)

  # Build the formula for the hierarchical model
  ranef = paste0("X", 2:n_trts, collapse = " + ")
  fixed = paste0(" + (1 + ", ranef ,"|id)")       # assumes id indexes participant
  model_formula = paste0("Y ~ ", ranef, fixed)

  # Generating posterior samples after data collected for everyone in period
  # Aggregate level analysis to get the population level effects
  posterior = rstanarm::stan_glmer(
    formula = model_formula,
    data = current_data, family = "gaussian",
    prior_intercept = rstanarm::normal(100, 10),
    prior = rstanarm::normal(0, 2.5),
    prior_aux = rstanarm::exponential(1, autoscale = TRUE),
    prior_covariance = rstanarm::decov(reg = 1, conc = 1, shape = 1, scale = 1),
    chains = chains, iter = iter, seed = 1, adapt_delta = adapt_delta,
    control = list(max_treedepth = max_treedepth),
    refresh = 0,
    QR = T
  )

  trial_posteriors[["burn"]] = as.data.frame(posterior)

  # Index to start the loop, since n_trts periods have been used in FRN phase
  adaptive_start = n_trts + 1 # FRN requires one period from every treatment

  # Calculate allocation probabilities based on Thompson Sampling
  # function for outputting allocation probabilities for each person
  # input: stan model, output: tibble of ids & allocation probabilites
  current_probs = allocate_probabilities(posterior, n_trts, optimize = optimize)
  current_probs$period = n_trts # By convention, just assign the last burn in period

  # Start collecting all of the allocation probabilities into one place
  trial_probs = current_probs

  # Run a treatment period for each person, updating the model after the
  # period has ended
  for (p in seq(adaptive_start, n_periods)) {

    print(paste0("Starting period ", p))
    trial_posteriors[[p]] = list()

    for (i in seq_len(n_subj)) {

      id_probs = current_probs %>% dplyr::filter(id == i)

      # Create treatment vector based on regime (1 where active, 0 else)
      next_trt = sample(1:n_trts, size = 1,
                        prob = id_probs %>% dplyr::select(starts_with("X")) %>%
                          unlist)
      x = matrix(0, nrow = 1, ncol = n_trts)
      x[1] = 1              # intercept
      x[next_trt] = 1       # setting next treatment
      design = x[rep(1, n_obvs), ]

      # Generate observations for subject based on selected treatment for cycle
      id_obvs = rep(x %*% betas[i, ], n_obvs)
      id_y = id_obvs + stats::rnorm(n_obvs, 0, y_sigma)

      if (n_obvs == 1) {
        id_data = tibble_row(id = i, Y = id_y, period = p)
        for (j in seq_len(n_trts)) {
          id_data[[paste0("X", j)]] = x[j]
        }
      } else {
        id_data = cbind(rep(i, n_obvs), design, id_y, rep(p, n_obvs)) %>% tibble::as_tibble()
        colnames(id_data) = c("id", paste0("X", 1:n_trts), "Y", "period")
      }

      current_data = dplyr::bind_rows(current_data, id_data)

      # Update the posterior distribution after observing these new values
      posterior = rstanarm::stan_glmer(
        formula = model_formula,
        data = current_data, family = "gaussian",
        prior_intercept = rstanarm::normal(100, 10),
        prior = rstanarm::normal(0, 2.5),
        prior_aux = rstanarm::exponential(1, autoscale = TRUE),
        prior_covariance = rstanarm::decov(reg = 1, conc = 1, shape = 1, scale = 1),
        chains = chains, iter = iter, seed = seed, adapt_delta = adapt_delta,
        control = list(max_treedepth = max_treedepth),
        refresh = 0,
        QR = T
      )

      trial_posteriors[[p]][[i]] = as.data.frame(posterior)

      # Add the probabilities used for the individual into the set
      # Need to do it now because it will possibly change with more data
      id_probs$period = p
      trial_probs = dplyr::bind_rows(trial_probs, id_probs)

      # Recalculate the allocation probabilities using the new data
      current_probs = allocate_probabilities(posterior, n_trts, optimize = optimize)

    } # end of subject loop

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
    data = current_data,
    posteriors = trial_posteriors
  )

  output

}
