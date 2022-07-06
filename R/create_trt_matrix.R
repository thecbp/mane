#' Internal function for constructing design matrix, used for Thompson Sampling
#'
#' @param n_trts Integer indicating number of treatments used in the trial
#'
#' @return Matrix indicating which treatment is currently active
#' @noRd
#'
#' @examples
#' n_trts = 3
#' T_mat = create_trt_matrix(n_trts)
create_trt_matrix = function(n_trts)  {

  mat = diag(n_trts)  # Diagonals represent treatment effects in use
  mat[,1] = 1         # intercept

  t(mat)
}
