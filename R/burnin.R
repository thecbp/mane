#' Title
#'
#' @param n_trts
#' @param n_burn_cycles
#' @param burn_obvs_per_period
#' @param betas
#' @param y_sigma
#' @param phi
#'
#' @return
#' @export
#'
#' @examples
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

  output = cbind(X, Y = Y) %>% tibble::as_tibble()
  output = output %>%
    dplyr::mutate(
      period = rep(1:n_trts, times = n_burn_cycles * burn_obvs_per_period)
    )

  output

}
