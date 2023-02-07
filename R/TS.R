TS = function(posterior, c, objective = "Maximum") {

  # Set the correct objective function to optimize for
  if (objective == "Maximum") { m = max }
  else { m = min }

  betas = as.data.frame(posterior) %>%
    dplyr::select(tidyselect::contains("b_"))

  n_trts = ncol(treatment_cols)

  # Design matrix representing outcome based on treatments
  X = diag(n_trts)
  X[,1] = 1
  colnames(X) = paste0("X", 1:n_trts)

  # Find which arm produced optimal outcome
  posterior_outcomes = X %*% t(betas)
  optimal_arms = apply(posterior_outcomes, 2, function(x) {which(x == m(x))})

  probs = c(table(optimal_arms) / nrow(betas))

  # Stabilize the probabilities
  stabilized_probs = probs^c / sum(probs^c)

  stabilized_probs

}
