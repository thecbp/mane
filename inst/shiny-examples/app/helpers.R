simulate = function(n_trts, n_burn_cycles, burn_obvs_per_period,
                    adaptive_obvs_per_period, max_duration,
                    betas, y_sigma, priors, n_chains, n_iter, lag,
                    stabilize = NULL,
                    objective = "Maximize",
                    adapt_delta = 0.999,
                    max_treedepth = 15,
                    seed = 1) {

  # Burn-in phase
  trial_data = burnin(n_trts, n_burn_cycles, burn_obvs_per_period, betas, y_sigma)

  # Build the formula for the model (X1 is reference treatment)
  f = paste0("Y~", paste0("X", 2:n_trts, collapse = "+"))

  # Setting priors for the model
  individual_priors = c(
    brms::set_prior(priors[["Intercept"]], class = "Intercept"),
    brms::set_prior(priors[["b"]], class = "b")
    # Using the default priors for sigma (half_t)
  )

  # Setting up controls for brms
  control = list(max_treedepth = max_treedepth,
                 adapt_delta = adapt_delta)

  burnin_periods = (n_trts * n_burn_cycles)

  # Generating posterior samples after data collected for everyone in period
  posteriors = brms::brm(formula = f,
                         data = trial_data,
                         family = gaussian(link = "identity"),
                         chains = n_chains,
                         iter = n_iter,
                         prior = individual_priors,
                         control = control,
                         refresh = 0)


  # Record all of the models in list, index is the period number
  trial_posteriors = list()
  trial_posteriors[[burnin_periods]] = posteriors

  # If stabilizing parameter not specified, use Thall & Walthen
  if (is.null(stabilize)) { c = burnin_periods / (2 * max_duration) }
  else { c = stabilize }

  # Calculate allocation probabilities via on Thompson Sampling
  current_probs = TS(posteriors, c, objective = objective)

  # Start collecting all of the allocation probabilities into one place
  trial_probs = c(current_probs, burnin_periods)

  # Index to start the loop, since n_trts periods have been used in FRN phase
  adaptive_start = burnin_periods + 1 # FRN requires one period from every treatment

  # Start of the adaptive phase of the trial
  for (per in adaptive_start:max_duration) {

    # Select the next treatment based on reallocated probabilities
    next_trt = sample(1:n_trts, size = 1, prob = current_probs)

    # Construct design matrix and "observe" the outcome under this treatment
    obv = matrix(0, nrow = 1, ncol = n_trts)
    obv[1] = 1        # intercept
    obv[next_trt] = 1 # setting next treatment
    X = obv[rep(1, adaptive_obvs_per_period), ]
    Y = (X %*% betas) + stats::rnorm(adaptive_obvs_per_period, 0, y_sigma)
    period_col = matrix(rep(per, adaptive_obvs_per_period), nrow = adaptive_obvs_per_period)
    new_data = cbind(X, Y, period_col)

    # Rename this new data to make sure that columns are aligned when combining
    colnames(new_data) = c(paste0("X", 1:n_trts), "Y", "period")

    trial_data = rbind(trial_data, new_data)

    # Recalculate the posteriors using the updated dataset
    # posteriors = brms::brm(formula = f,
    #                        data = trial_data,
    #                        family = gaussian(link = "identity"),
    #                        chains = n_chains,
    #                        iter = n_iter,
    #                        prior = individual_priors,
    #                        control = control,
    #                        refresh = 0)
    posteriors = update(posteriors, newdata = trial_data)

    # Record and update the posteriors used in the loop
    trial_posteriors[[per]] = posteriors

    # If stabilize not initialized, use Thall & Walthen's optimal parameter
    if (is.null(stabilize)) { c = per / (2 * max_duration) }
    else { c = stabilize }

    # Readjust reallocation probabilities based on new updated posterior
    current_probs = TS(posteriors, c, objective = objective)

    # Record the new probabilities & period
    trial_probs = rbind(trial_probs, c(current_probs, per))

  } # end of period loop

  # Post simulation processing
  trial_probs = tibble::as_tibble(trial_probs)
  names(trial_probs) = c(paste0("X", 1:n_trts), "period")



  # Output of simulation
  out = list()
  out[["trial_params"]] = list(
    n_trts = n_trts,
    n_burn_cycles = n_burn_cycles,
    burn_obvs_per_period = burn_obvs_per_period,
    adaptive_obvs_per_period = adaptive_obvs_per_period,
    max_duration = max_duration,
    betas = betas,
    y_sigma = y_sigma,
    priors = priors
  )
  out[["stan_params"]] = list(
    n_chains = n_chains,
    n_iter = n_iter,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )
  out[["data"]] = trial_data
  out[["posteriors"]] = trial_posteriors
  out[["allocation_probs"]] = trial_probs

  out



}
