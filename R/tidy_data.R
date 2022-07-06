#' Function for converting Stan-friendly data to be more tidyverse friendly
#'
#' @param standata List containing the data that was
#' @param n_obvs not needed?
#'
#' @return
#' @internal
#'
#' @examples
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
