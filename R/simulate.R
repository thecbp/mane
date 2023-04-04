#' Function for simulating a single Platform-of-1 design
#'
#' @param n_trts (integer): The number of treatments. The value of \code{n_trts} should be a positive integer.
#' @param n_burn_cycles (integer): The number of burn-in cycles. This parameter defines how many times the dataset will be replicated for the burn-in phase. The value of \code{n_burn_cycles} should be a positive integer.
#' @param burn_obvs_per_period (integer): The number of observations per period for each treatment during the burn-in period. The value of \code{burn_obvs_per_period} should be a positive integer.
#' @param adaptive_obvs_per_period (integer): The number of observations per period for each treatment during the adaptive phase of the trial. The value of \code{adaptive_obvs_per_period} should be a positive integer.
#' @param max_duration (integer): The maximum duration of the trial, expressed as the number of periods. The value of \code{max_duration} should be a positive integer.
#' @param betas (numeric vector): A vector of regression coefficients representing the effects of the treatments. The length of this vector should be equal to \code{n_trts}.
#' @param y_sigma (numeric): The standard deviation of the error term in the outcome variable. This parameter represents the noise in the outcome variable and should be a positive number.
#' @param priors  (list): A list of priors to be used in the model. See the example to see how to format this list.
#' @param n_chains (integer): The number of chains to use for MCMC. The value of n_chains should be a positive integer.
#' @param n_iter (integer): The number of iterations for MCMC. The value of n_iter should be a positive integer.
#' @param phi (numeric): The autoregressive (AR) coefficient for the outcome variable. If phi is 0, there is no serial correlation in the outcome variable. Otherwise, the outcome variable is generated with an AR(1) process using the specified phi value.
#' @param stabilize (numeric, optional): The stabilizing parameter for the trial. If not provided, the function uses Thall & Walthen's optimal parameter. The value of stabilize should be a positive number or NULL.
#' @param objective (string): The objective of the trial, either "Maximize" or "Minimize". Determines the goal of the Thompson Sampling algorithm.
#' @param adapt_delta (numeric): The adaptation parameter for the Hamiltonian Monte Carlo algorithm used in the brms package. The value of \code{adapt_delta} should be between 0 and 1.
#' @param max_treedepth  (integer): The maximum tree depth for the Hamiltonian Monte Carlo algorithm used in the brms package. The value of \code{max_treedepth} should be a positive integer.
#'
#' @return A list containing the following elements:
#' - \code{trial_params}: A list of the input trial parameters.
#' - \code{stan_params}: A list of the input Stan parameters (used in the brms package).
#' - \code{data}: A data frame containing the simulated trial data.
#' - \code{burnin_posteriors}: A brmsfit that contains the posteriors after the burn-in phase
#'
#' @export
#'
#' @examples
#' # Simulate a single Platform-of-1 trial
#' sims = simulate(n_trts = 3,
#'                 n_burn_cycles = 1,
#'                 burn_obvs_per_period = 7,
#'                 adaptive_obvs_per_period = 7,
#'                 max_duration = 12,
#'                 betas = c(130, -10, -5),
#'                 y_sigma = 10,
#'                 priors = list("Intercept" = "normal(0,100)",
#'                               "b" = "normal(0,100)"),
#'                 phi = 0,
#'                 objective = "Minimize")
simulate = function(n_trts, n_burn_cycles, burn_obvs_per_period,
                    adaptive_obvs_per_period, max_duration,
                    betas, y_sigma, priors, n_chains, n_iter, phi,
                    stabilize = NULL,
                    objective = "Maximize",
                    adapt_delta = 0.999,
                    max_treedepth = 15) {



  # Burn-in phase
  if (phi == 0) {
    trial_data = burnin(n_trts, n_burn_cycles, burn_obvs_per_period, betas, y_sigma)
  } else {
    trial_data = burnin_corr(n_trts, n_burn_cycles, burn_obvs_per_period, betas, y_sigma, phi = phi)
  }

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

  n_cores = parallel::detectCores()

  # Generating posterior samples after data collected for everyone in period
  posteriors = brms::brm(formula = f,
                         data = trial_data,
                         family = gaussian(link = "identity"),
                         chains = n_chains,
                         iter = n_iter,
                         prior = individual_priors,
                         control = control,
                         cores = n_cores,
                         refresh = 0)

  burnin_posteriors = posteriors

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
    X = obv[rep(1, adaptive_obvs_per_period), ] %>% matrix(nrow = adaptive_obvs_per_period)

    # Generate outcome based on if we need serial correlation or not
    if (phi == 0) {
      Y = (X %*% betas) + stats::rnorm(adaptive_obvs_per_period, 0, y_sigma)
    } else {
      Y_ar = arima.sim(model = list(ar = phi), sd = y_sigma, n = adaptive_obvs_per_period)
      Y = Y_ar + (X %*% betas)
    }

    period_col = matrix(rep(per, adaptive_obvs_per_period), nrow = adaptive_obvs_per_period)

    new_data = cbind(X, Y, period_col)

    # Rename this new data to make sure that columns are aligned when combining
    colnames(new_data) = c(paste0("X", 1:n_trts), "Y", "period")

    trial_data = rbind(trial_data, new_data)

    # Recalculate the posteriors using the updated dataset
    posteriors = update(posteriors, newdata = trial_data)

    # If stabilize not initialized, use Thall & Walthen's optimal parameter
    if (is.null(stabilize)) { c = per / (2 * max_duration) }
    else { c = stabilize }

    # Readjust reallocation probabilities based on new updated posterior
    current_probs = TS(posteriors, c, objective = objective)

    # Record the new probabilities & period
    trial_probs = rbind(trial_probs, c(current_probs, per))

  } # end of adaptive phase loop

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
    phi = phi,
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
  out[["burnin_posteriors"]] = burnin_posteriors
  out[["allocation_probs"]] = trial_probs

  out

}
