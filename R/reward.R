#' Function for translating posterior parameters into reward
#'
#' @param params vector of doubles representing parameters of the model
#' @param data matrix representing the treatment
#'
#' @return
#' @export
#'
#' @examples
reward = function(params, data) {
  params %*% data
}
