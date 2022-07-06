#' Function for stabilizing allocation probabilities
#'
#' @param probs
#' @param c
#'
#' @return
#' @internal
#'
#' @examples
stabilize_probs = function(probs, c) {
  # Stabilize the allocation probability according to Thall & Wathen

  probs_c = probs^c
  probs_c / sum(probs_c)

}
