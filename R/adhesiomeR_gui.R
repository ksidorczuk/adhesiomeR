#' Run adhesiomeR GUI
#' 
#' Launches a graphical user interface to adhesiomeR in a form of
#' a shiny application. The application may be used to perform adhesiomeR
#' analysis on up to 100 genomes.
#' @return No return value, called for its side effects.
#' @importFrom shiny runApp
adhesiomeR_gui <- function() {
  runApp(system.file("adhesiomeR", package = "adhesiomeR"))
}