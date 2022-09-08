#' Reallocate the treatment randomization probabilities via Thompson Sampling
#'
#' @param posterior stanreg object that contains the posterior samples
#' @param n_trts integer indicating number of treatments currently in trial
#' @param optimize character indicating if the maximum or minimum reward should be optimized for
#' @param level character indicating if the model is individual or aggregate-level
#'
#' @return Dataframe containing the treatment allocation probabilities by id
#' @export
#'
allocate_probabilities = function(posterior,
                                  n_trts,
                                  optimize = "max",
                                  level = "aggregate") {

  # Place posterior samples in a data frame
  samples = as.data.frame(posterior)
  S = nrow(samples)

  # construct prob_df so that it initializes an empty matrix of zeroes for each trt

  if (optimize == "max") { m = max}
  else { m = min }

  if (level == "aggregate") {

    # Back-calculating number of subjects from posterior dataframe
    n_subj = (ncol(samples) - (n_trts) - ((n_trts * (n_trts+1) * 0.5)) - 1) * (1/n_trts)
    prob_df = matrix(0, nrow = n_subj, ncol = n_trts)

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
      find_max = function(x) {which(x == m(x))}
      optimal_arms = apply(posterior_outcomes, 2, find_max)

      probs = c(table(optimal_arms)) / S

      prob_df[i, as.numeric(names(probs))] = probs
      colnames(prob_df) = paste0("X", 1:n_trts)

    }

  } else if (level == "individual") {

    n_subj = (ncol(samples) - ((n_trts * (n_trts+1) * 0.5)) - 1) * (1/n_trts)
    prob_df = matrix(0, nrow = n_subj, ncol = n_trts)

    # Iterate over all of the subjects
    for (i in seq_len(n_subj)) {

      # Design matrix representing outcome based on treatments
      trt = diag(n_trts)
      trt[, 1] = 1
      colnames(trt) = paste0("X", 1:n_trts)

      # Calculate subject-specific betas
      current_id = paste0("id:",i, "]")
      betas = samples %>% dplyr::select(ends_with(current_id))

      # Find which arm produced best treatment effects
      posterior_outcomes = trt %*% t(betas)
      find_max = function(x) {which(x == m(x))}
      optimal_arms = apply(posterior_outcomes, 2, find_max)

      probs = c(table(optimal_arms)) / S

      prob_df[i, as.numeric(names(probs))] = probs
      colnames(prob_df) = paste0("X", 1:n_trts)

    }


  }

  output = prob_df %>%
    tibble::as_tibble() %>%
    mutate(id = 1:n_subj) %>%
    select(id, dplyr::everything())

  output

  # Potential TODO: place the stabilization here as well
  # Note: potentially wasteful here because this is done for everyone each time
}
