#' Title
#'
#' @return
#' @export
#'
#' @examples
runPlanner = function() {
  appDir = system.file("shiny-examples", "app", package = "mane")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mane`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
