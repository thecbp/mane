#' Internal function for stabilizing allocation probabilities in Thompson Sampling
#'
#' Stabilize the allocation probability according to Thall & Wathen
#'
#' @param probs Vector of doubles indicating treatment allocation probabilities
#' @param c Tuning parameter for stabilizing the probabilities
#'
#' @return Vector of doubles indicating treatment allocation probabilities
#' @noRd
#'
#' @examples
#' probs = runif(3)
#' stabilized_probs = stabilize_probs(probs, 0.75)
stabilize_probs = function(probs, c) {

  probs_c = probs^c
  probs_c / sum(probs_c)

}
