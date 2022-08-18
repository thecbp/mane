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

  # Calculate outcome based on each of the treatments
  trt = diag(n_trts)
  trt[, 1] = 1
  colnames(trt) = paste0("X", 1:n_trts)

  # first n_subj elements come from first person
  Y = (trt %*% t(betas)) %>% c()

  # Replicate treatment matrix to number of subjects
  trts_by_subj = trt[rep(seq_len(nrow(trt)), n_subj),]

  # Append the outcome
  trts_by_subj = cbind(id = rep(1:n_subj, each = n_trts),
                       period = rep(1:n_subj, times = n_trts),
                       trts_by_subj,
                       Y = Y)

  # Duplicate the matrix to match number of observations for each person
  full_trt_mat = trts_by_subj[rep(seq_len(nrow(trts_by_subj)), n_obvs),] %>%
    as_tibble()

  # Add subject level noise
  full_trt_mat$Y = full_trt_mat$Y + rnorm(nrow(full_trt_mat), 0, y_sigma)

  full_trt_mat
}
