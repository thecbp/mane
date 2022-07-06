#' Simulate a fixed randomization design for an N-of-1 trial
#'
#' @param n_subj Number of subjects to include in trial
#' @param n_trts Number of treatments to be considered
#' @param n_periods Number of treatments periods
#' @param n_obvs Number of observations made during a treatment period
#' @param betas Treatment effects for each of the n_trts treatments
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
#' current_data = generate_FRN_data(n_subj = 2, n_trts = 3, n_cycles = 1,
#'                                  n_obvs = 5, betas, y_sigma = 2)
generate_FRN_data = function(n_subj,
                             n_trts,
                             n_periods,
                             n_obvs,
                             betas,
                             y_sigma) {

  # Setting up parameters
  N = n_subj * n_periods * n_trts * n_obvs # total sample size of overall trial
  y = id = period = array(0, N)            # initializing data structures
  X = array(0, N * n_trts) %>%
    matrix(nrow = N, ncol = n_trts)        # Set up data matrix for treatments
  iter = 1                                 # global counter

  # Simulate data from fixed randomization scheme
  # For each treatment cycle...
  for (p in 1:n_periods) {

    # ... for each subject...
    for (i in 1:n_subj){

      # ... randomize treatment order for the cycle...
      regime = sample(1:n_trts, size = n_trts, replace = F)

      # ... for each of their treatment periods...
      for (trt in 1:n_trts) {

        # ... measure their outcome [n_obvs] times
        for (z in 1:n_obvs) {

          # Create treatment vector based on regime (1 where active, 0 else)
          x = array(0, n_trts)
          x[1] = 1              # intercept
          x[regime[trt]] = 1    # setting treatment, if not control

          X[iter,] = x                                  # treatment vector
          y[iter] = x %*% betas[i,] + stats::rnorm(1, 0, y_sigma)  # outcome
          id[iter] = i                                  # subject ID
          period[iter] = trt + n_trts * (p - 1)         # period number

          iter = iter + 1
        }
      }
    }
  }

  # Make data usable for Stan
  list(J = n_subj,
       K = n_trts,
       N = N,
       X = X,
       y = y,
       id = id,
       period = period)
}
