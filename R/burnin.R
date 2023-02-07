burnin = function(n_trts, n_burn_cycles, burn_obvs_per_period, betas, y_sigma) {

  # Matrix of single person getting all of the treatments for n_obvs each
  X = diag(n_trts)
  X[, 1] = 1
  colnames(X) = paste0("X", 1:n_trts)
  X = X[rep(seq_len(nrow(X)), each = burn_obvs_per_period),]
  X = X[rep(seq_len(nrow(X)), times = n_burn_cycles),]

  # Calculate outcome with linear model
  # Also adding within-subject noise
  Y = ((X %*% t(betas)) %>% c()) + stats::rnorm(nrow(X), 0, y_sigma)

  output = cbind(X, Y = Y) %>% tibble::as_tibble()
  output = output %>%
    dplyr::mutate( period = rep(1:n_trts, each = n_obvs) ) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate( time = 1:n_obvs) # Mark the time of observation

  output

}
