#' Reallocate the treatment randomization probabilities via Thompson Sampling
#'
#' @param posterior stanreg object that contains the posterior samples
#' @param c float between 0 and 1 indicating how much to stablize probabilies
#' @param objective character indicating whether to optimize min or max

#'
#' @return Dataframe containing the treatment allocation probabilities by id
#' @export
#'
allocate_probabilities = function(posterior, c, objective = "max") {

  # Set the correct function to maximize on
  if (objective == "max") { m = max }
  else { m = min }

  output = posterior %>%
    dplyr::transmute(
      id = id,
      probs = purrr::map(model, function(mod) {

        betas = as.data.frame(mod) %>% dplyr::select(tidyselect::contains("b_"))
        n_trts = ncol(betas)

        # Design matrix representing outcome based on treatments
        trt = diag(n_trts)
        trt[, 1] = 1
        colnames(trt) = paste0("X", 1:n_trts)

        # Find which arm produced best treatment effects
        posterior_outcomes = trt %*% t(betas)
        find_optimal = function(x) {which(x == m(x))}
        optimal_arms = apply(posterior_outcomes, 2, find_optimal)

        probs = c(table(optimal_arms) / nrow(betas))

        # Stabilize the probabilities
        stab_probs = probs^c / sum(probs^c)

        prob_df = stab_probs %>%
          tibble::as_tibble() %>%
          dplyr::mutate(
            trt = names(stab_probs)
          ) %>%
          tidyr::pivot_wider(names_from = trt, values_from = value)

        prob_df
      })
    ) %>%
    tidyr::unnest(probs)

  # Adjust column names
  n_trts = ncol(output) - 1
  colnames(output) = c("id", paste0("X", 1:n_trts))

  # Replace NA with 0
  replace_list = list()
  for (col in paste0("X", 1:n_trts)) {
    replace_list[[col]] = 0
  }

  output = output %>% replace_na(replace_list)

  output
}
