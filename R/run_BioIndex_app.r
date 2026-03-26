#' Run the BioIndex Shiny application
#'
#' Launches the embedded Shiny application included in the package.
#'
#' @export
run_BioIndex_app <- function() {
    app_dir <- system.file("shiny/BioIndexApp", package = "BioIndex")

    if (app_dir == "") {
        stop("Could not find the BioIndex Shiny app directory. Please reinstall the package.")
    }

    shiny::runApp(app_dir, display.mode = "normal")
}
