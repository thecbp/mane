#' Function for simulating the burn-in phase of the Platform-of-1 design
#'
#' @param n_trts (integer): The number of treatments. The value of n_trts should be a positive integer.
#' @param n_burn_cycles (integer): The number of burn-in cycles. This parameter defines how many times the dataset will be replicated for the burn-in phase. The value of n_burn_cycles should be a positive integer.
#' @param burn_obvs_per_period  (integer): The number of observations per period for each treatment during the burn-in period. The value of burn_obvs_per_period should be a positive integer.
#' @param betas (numeric vector): A vector of regression coefficients representing the effects of the treatments. The length of this vector should be equal to n_trts.
#' @param y_sigma (numeric): The standard deviation of the error term in the outcome variable. This parameter represents the noise in the outcome variable and should be a positive number.
#' @param phi (numeric): The autoregressive (AR) coefficient for the outcome variable. If phi is 0, there is no serial correlation in the outcome variable. Otherwise, the outcome variable is generated with an AR(1) process using the specified phi value.
#'
#' @return A tibble containing treatment indicators for all treatments and an outcome
#' @export
burnin = function(n_trts, n_burn_cycles, burn_obvs_per_period, betas, y_sigma, phi) {

  # Matrix of single person getting all of the treatments for n_obvs each
  X = diag(n_trts)
  X[, 1] = 1
  colnames(X) = paste0("X", 1:n_trts)
  X = X[rep(seq_len(nrow(X)), each = burn_obvs_per_period),]
  X = X[rep(seq_len(nrow(X)), times = n_burn_cycles),]

  # Generate outcome based on if we need serial correlation or not
  if (phi == 0) {
    Y = ((X %*% betas) %>% c()) + stats::rnorm(nrow(X), 0, y_sigma)
  } else {
    Y_ar = arima.sim(model = list(ar = phi), sd = y_sigma,
                     n = n_trts * n_burn_cycles * burn_obvs_per_period)
    Y = Y_ar + X %*% betas
  }

  output = cbind(X, Y = Y) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      period = rep(1:n_trts, times = n_burn_cycles * burn_obvs_per_period)
    )

  output

}
