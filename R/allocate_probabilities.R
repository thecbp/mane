#' Reallocate the treatment randomization probabilities via Thompson Sampling
#'
#' @param posterior stanreg object that contains the posterior samples
#' @
#'
#' @return
#' @export
#'
#' @examples
allocate_probabilities = function(posterior, n_trts) {

  # Place posterior samples in a data frame
  samples = as.data.frame(posterior)
  S = nrow(samples)

  # Back-calculating number of subjects from posterior dataframe
  n_subj = (ncol(samples) - (n_trts) - ((n_trts * (n_trts+1) * 0.5)) - 1) * (1/n_trts)

  prob_df = tibble()

  # Iterate over all of the subjects
  for (i in seq_len(n_subj)) {

    # Design matrix representing outcome based on treatments
    trt = diag(n_trts)
    trt[, 1] = 1
    colnames(trt) = paste0("X", 1:n_trts)

    # Calculate subject-specific betas
    subj_ranef_index = (1:n_trts) + (n_trts * i)
    betas = samples[,1:n_trts] + samples[,subj_ranef_index]

    # Find which arm produced best treatment effects
    posterior_outcomes = trt %*% t(betas)
    find_max = function(x) {which(x == max(x))}
    optimal_arms = apply(posterior_outcomes, 2, find_max)

    probs = c(table(optimal_arms)) / S
    prob_df = bind_rows(prob_df, c(id = i, probs))

  }

  fill_list = list()
  for (t in seq_len(n_trts)) fill_list[[as.character(t)]] = 0
  prob_df %>% replace_na(replace = fill_list)


  # Potential TODO: place the stabilization here as well
  # Note: potentially wasteful here because this is done for everyone each time
}
