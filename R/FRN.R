#' Simulate a fixed randomization design for a multi-arm N-of-1 trial
#'
#' @param n_subj Integer indicating number of subjects in trial
#' @param n_trts Integer indicating number of treatments in trial
#' @param n_periods
#' @param n_obvs
#' @param betas
#' @param y_sigma
#' @param chains
#' @param warmup
#' @param iter
#' @param adapt_delta
#' @param max_treedepth
#'
#' @return
#' @export
#'
#' @examples
FRN = function(n_subj, n_trts, n_periods, n_obvs, betas, y_sigma,
               chains, warmup, iter, adapt_delta = 0.99, max_treedepth = 15) {

  # Fixed randomization N-of-1 posterior sampling

  # n_subj:
  # n_trts: number of treatments
  # n_cycles: number of cycles
  # n_obvs: number of observations per period
  # betas: matrix containing regression coefficients for outcome by subject
  # y_sigma: within-person standard deviation (noise in outcome model)
  # chains: number of chains to use in MCMC
  # warmup: number of observations to use to warmup the MCMC
  # iter: how long the chains should go total

  # Entire trial is run on fixed randomization scheme
  current_data = generate_FRN_data(n_subj,
                                   n_trts,
                                   n_periods,
                                   n_obvs,
                                   betas,
                                   y_sigma)

  # Make the data easier to visualize
  trial_data = tidy_data(current_data, n_obvs = n_obvs)

  # Generating posterior samples after data collected for everyone in period
  stan_model = mcmc_agg(current_data, chains = chains, warmup = warmup, iter = iter,
                        adapt_delta = adapt_delta, max_treedepth = max_treedepth)

  list(
    n_subj = n_subj,
    n_periods = n_periods,
    n_obvs = n_obvs,
    betas = betas,
    y_sigma = y_sigma,
    data = trial_data,
    stan_model = stan_model
  )


}

