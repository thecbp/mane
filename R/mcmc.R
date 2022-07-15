#' Wrapper for running an rstan model and getting posterior samples
#'
#' @param stanfile Character vector indicating path of Stan model (.stan)
#' @param data List of data from the trial
#' @param chains Integer indicating number of chains to use in MCMC
#' @param warmup Integer indicating how long the warmup period for MCMC should be
#' @param iter Integer indicating the total number of posterior samples to generate
#' @param adapt_delta Integer indicating the step size for MCMC
#' @param max_treedepth Integer indicating the maximum tree depth for MCMC
#'
#' @return List containing posterior samples of treatment effect & trial parameters
#' @noRd
#'
#' @examples
#' stanfile = 'aggregated.stan'
#' data = generate_FRN_data(n_subj = 2, n_trts = 3, n_periods = 3,
#'                          n_obvs = 5, betas, y_sigma = 2)
#' mcmc = mcmc(stanfile, data)
mcmc = function(stanfile, data, chains = 1, warmup = 1000, iter = 3000,
                adapt_delta = 0.99, max_treedepth = 15) {

  h_out = rstan::stan(stanfile,
                      data = data,
                      chains = chains,         # number of Markov chains
                      warmup = warmup,         # number of warmup iterations per chain
                      iter = iter,             # total number of iterations per chain
                      refresh = 0,
                      verbose = TRUE,
                      control = list(adapt_delta = adapt_delta,
                                     max_treedepth = max_treedepth))

  h_out
}
