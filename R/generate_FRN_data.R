#' Simulate a fixed randomization design for an N-of-1 trial
#'
#' @param n_subj Number of subjects to include in trial
#' @param n_trts Number of treatments to be considered
#' @param n_periods Number of treatments periods
#' @param n_obvs Number of observations made during a treatment period
#' @param betas Matrix containing treatment effects for each individual (rows) and treatment (columns)
#' @param y_sigma Within-individual noise of continuous outcome
#'
#' @return A list with design parameters and observed treatment and outcomes
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # Generate array of treatment effects
#' # Rows are subjects, columns are treatments
#' betas = array(1:6, dim = c(2, 3))
#' current_data = generate_FRN_data(n_subj = 2, n_trts = 3, n_periods = 3,
#'                                  n_obvs = 5, betas, y_sigma = 2)
generate_FRN_data = function(n_subj,
                             n_trts,
                             n_periods,
                             n_obvs,
                             betas,
                             y_sigma) {

  # Matrix of single person getting all of the treatments for n_obvs each
  trt = diag(n_trts)
  trt[, 1] = 1
  colnames(trt) = paste0("X", 1:n_trts)
  trt = trt[rep(seq_len(nrow(trt)), each = n_obvs),] # single design matrix
  X = trt[rep(seq_len(nrow(trt)), each = n_subj),] # matrices for each person

  # Calculate outcome with given betas (column is person, row is outcome for first obvs)
  # Also adding within-subject noise
  Y = ((trt %*% t(betas)) %>% c()) + stats::rnorm(nrow(X), 0, y_sigma)

  id = rep(1:n_subj, times = n_obvs * n_trts)

  output = cbind(id = id, X, Y = Y) %>% tibble::as_tibble()

  output
}
