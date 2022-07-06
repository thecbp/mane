#' Wrapper for running rstan
#'
#' @param data List of data from the trial
#' @param chains
#' @param warmup
#' @param iter
#' @param adapt_delta
#' @param max_treedepth
#'
#' @return
#' @internal
#'
#' @examples
mcmc = function(file, data, chains = 1, warmup = 1000, iter = 3000,
                adapt_delta = 0.99, max_treedepth = 15) {

  h_out = stan(file = file, # file for aggregated analysis
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
