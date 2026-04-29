#' launches the embedded Shiny application included in the package.
#'
#' @author Walter Zupa \email{zupa@fondazionecoispa.org}
#'
#' @importFrom shiny runApp
#' @return \code{No return value, called for side effects} (launches the Shiny application).
#' @export
run_BioIndex_app <- function() {
    app_dir <- system.file("shiny/BioIndexApp", package = "BioIndex")

    if (app_dir == "") {
        stop("Could not find the BioIndex Shiny app directory. Please reinstall the package.")
    }

    shiny::runApp(app_dir, display.mode = "normal")
}
