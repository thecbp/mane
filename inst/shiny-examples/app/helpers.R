simulate = function(n_trts, n_burn_cycles, burn_obvs_per_period,
                    adaptive_obvs_per_period, max_duration,
                    betas, y_sigma, priors, n_chains, n_iter, lag,
                    stabilize = NULL,
                    objective = "maximize",
                    adapt_delta = 0.999,
                    max_treedepth = 15,
                    seed = 1) {

  set.seed(seed)

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

  # DELETE LATER: burnin() outputs a tibble: X1, X2, X3, Y, time

  # Generating posterior samples after data collected for everyone in period
  posteriors = brms::brm(data = trial_data, formula = f,
                         family = gaussian(link = "identity"),
                         chains = n_chains,
                         iter = n_iter,
                         prior = individual_priors,
                         control = control,
                         refresh = 0)


  # Record all of the models in list, index is the period number
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
    Y = X %*% betas + stats::rnorm(n_obvs, 0, y_sigma)

    new_data = cbind(X, Y)
    new_data$period = per

    colnames(new_data) = c(paste0("X", 1:n_trts), "Y", "period")

    trial_data = rbind(current_data, new_data)
    trial_data$time = 1:nrow(current_data)

    # Recalculate the posteriors using the updated dataset
    posteriors = brms::brm(data = trial_data, formula = f,
                           family = gaussian(link = "identity"),
                           chains = n_chains,
                           iter = n_iter,
                           prior = individual_priors,
                           control = control,
                           refresh = 0)

    # Record and update the posteriors used in the loop
    trial_posteriors[[per]] = posteriors

    # If stabilize not initialized, use Thall & Walthen's optimal parameter
    if (is.null(stabilize)) { c = per / (2 * n_periods) }
    else { c = stabilize }

    # Readjust reallocation probabilities based on new updated posterior
    current_probs = TS(posteriors, c, objective = objective)

    # Record the new probabilities & period
    trial_probs = rbind(trial_probs,
                        c(current_probs, per))

  } # end of period loop

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
