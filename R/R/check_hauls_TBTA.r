#' Check hauls TB TA (RoME)
#' @description
#' Check if all the hauls in TB are in TA.
#'
#' @param DataTA MEDITS TA table
#' @param DataTB MEDITS TB table
#' @param year reference year for the analysis
#' @param wd working directory
#' @param suffix name of the log file
#' @note
#' This function is an internal routine based on \strong{RoME version 0.2.3}.
#' It is provided within \code{BioIndex} to ensure the package remains
#' functional and self-sufficient for data validation.
#' @export
check_hauls_TBTA <- function(DataTA, DataTB, year, wd=NA, suffix=NA) {

  if (is.na(wd)) {
    wd <- tempdir()
  }

  if (!file.exists(file.path(wd, "Logfiles"))) {
    dir.create(file.path(wd, "Logfiles"), recursive = TRUE, showWarnings = FALSE)
  }
  
  if (missing(suffix) || is.na(suffix)) {
    suffix <- paste0(Sys.Date(), format(Sys.time(), "_time_h%Hm%Ms%OS0"))
  }
  
  Errors <- file.path(wd, "Logfiles", paste0("Logfile_", suffix, ".dat"))
  if (!file.exists(Errors)) {
    file.create(Errors)
  }

  if (length(year) != 1 || is.na(year)) {
    stop("Argument 'year' must be a single non-NA value")
  }
  
  DataTA <- DataTA[DataTA$YEAR == year, ]
  DataTB <- DataTB[DataTB$YEAR == year, ]

  if (nrow(DataTB) == 0) {
    write(paste("\n----------- check presence in TA of TB hauls - ", year), file = Errors, append = TRUE)
    write("No data found for this year in TB", file = Errors, append = TRUE)
    return(TRUE)
  }

  write(paste("\n----------- check presence in TA of TB hauls - ", year), file = Errors, append = TRUE)

  numberError <- 0
  unique_hauls_TB <- unique(DataTB$HAUL_NUMBER)
  
  for (h in unique_hauls_TB) {
    if (!(h %in% DataTA$HAUL_NUMBER)) {
      write(paste("No haul", h, "in TA"), file = Errors, append = TRUE)
      numberError <- numberError + 1
    }
  }

  if (numberError == 0) {
    write("No error occurred", file = Errors, append = TRUE)
    return(TRUE)
  } else {
    return(FALSE)
  }
}
