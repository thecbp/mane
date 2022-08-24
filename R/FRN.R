#' Simulate a fixed randomization design for a multi-arm N-of-1 trial
#'
#'  `FRN()` simulates a fixed randomization scheme for an N-of-1 trial. The
#'  treatments are randomly ordered and given to the subject. This process is
#'  repeated for some number of cycles. At the end of the trial simulation, an
#'  analysis on the population level treatment effect is performed.
#'
#' @param n_subj Integer indicating number of subjects in trial
#' @param n_trts Integer indicating number of treatments in trial
#' @param n_periods Integer indicating number of treatment periods in trial
#' @param n_obvs Integer indicating how many observations should be recorded from a particular treatment regime
#' @param betas Matrix containing treatment effects for each individual (rows) and treatment (columns)
#' @param y_sigma Double indicating amount of within-person variance
#' @param chains Integer indicating number of chains to use in MCMC
#' @param warmup Integer indicating how long the warmup period for MCMC should be
#' @param iter Integer indicating the total number of posterior samples to generate
#' @param adapt_delta Integer indicating the step size for MCMC
#' @param max_treedepth Integer indicating the maximum tree depth for MCMC
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
#'1
#' frn = FRN(n_subj = 2, n_cycles = 3, n_obvs = 5, n_trts = 3, y_sigma = 2,
#'           satnfile = stanfile, betas = betas, chains = 1, warmup = 1000,
#'           iter = 3000)
#' }
FRN = function(n_subj, n_trts, n_periods, n_obvs, betas, y_sigma,
               chains, warmup, iter, adapt_delta = 0.999, max_treedepth = 17) {



  # Entire trial is run on fixed randomization scheme
  current_data = generate_FRN_data(n_subj,
                                   n_trts,
                                   n_periods,
                                   n_obvs,
                                   betas,
                                   y_sigma)


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
    QR = T
  )

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
      warmup = warmup,
      iter = iter,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    ),
    data = current_data,
    posterior = posterior
  )

  output

}

