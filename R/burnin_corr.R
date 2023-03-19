burnin_corr = function(n_trts, n_burn_cycles, burn_obvs_per_period, betas, y_sigma, phi) {

  # Matrix of single person getting all of the treatments for n_obvs each
  X = diag(n_trts)
  X[, 1] = 1
  colnames(X) = paste0("X", 1:n_trts)
  X = X[rep(seq_len(nrow(X)), each = burn_obvs_per_period),]
  X = X[rep(seq_len(nrow(X)), times = n_burn_cycles),]

  # Calculate outcome according to AR(p)
  Y_ar = arima.sim(model = list(ar = phi), sd = y_sigma,
                n = n_trts * n_burn_cycles * burn_obvs_per_period)

  Y = Y_ar + X %*% betas

  output = cbind(X, Y = Y) %>% tibble::as_tibble()
  colnames(output) = c(colnames(X), "Y")
  output = output %>%
    dplyr::mutate(
      period = rep(1:n_trts, times = n_burn_cycles * burn_obvs_per_period)
    )

  output

}
