#' Launch Shiny App for MultiCoex
#'
#' A function that launches the Shiny app for MultiCoex, an interactive
#' application to run the package. The code has been placed in
#' \code{./inst/shiny-scripts}.
#'
#' @return No return value but will open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' MultiCoex::runMultiCoex()
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. \href{https://shiny.rstudio.com/tutorial/}{Link}
#' Anjali Silva
#'
#' @export
#' @importFrom shiny runApp

runTestingPackage <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "MultiCoex")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")

  return(actionShiny)
}
# [END]
