#' Function for converting Stan-friendly data to be more tidyverse friendly
#'
#' @param standata List containing the data that was
#' @param n_obvs not needed?
#'
#' @return
#' @noRd
#'
#' @example
#' data = generate_FRN_data(n_subj = 2, n_trts = 3, n_periods = 3,
#'                          n_obvs = 5, betas, y_sigma = 2)
#' tidied = tidy_data(data)
tidy_data = function(standata) {

  data = standata$X %>% as.data.frame() %>% as_tibble()
  N = standata$N
  J = standata$J
  n_trts = standata$K
  colnames(data) = paste0("X", 1:n_trts)

  data = bind_cols(data,
                   id = standata$id,
                   y = standata$y,
                   period = standata$period) %>% unique

  return(data)

}
