#' Title
#'
#' @param n_trts
#'
#' @return
#' @internal
#'
#' @examples
create_trt_matrix = function(n_trts)  {
  # structure for quickly figuring out optimal arm, given betas

  mat = diag(n_trts)  # Diagonals represent treatment effects in use
  mat[,1] = 1         # intercept

  t(mat)
}
