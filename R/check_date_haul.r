#' Check date haul (RoME)
#' @description
#' Check if in TB, TC or TE the date by haul is the same of the one reported in TA.
#'
#' @param DataTA MEDITS TA table
#' @param Data MEDITS TB, TC or TE table
#' @param year reference year for the analysis
#' @param wd working directory
#' @param suffix name of the log file
#' @note
#' This function is an internal routine based on \strong{RoME version 0.2.3}.
#' It is provided within \code{BioIndex} to ensure the package remains
#' functional and self-sufficient for data validation.
#' @return \code{TRUE} if the validation passes, \code{FALSE} otherwise.
#' @export
check_date_haul <- function(DataTA, Data, year, wd=NA, suffix=NA) {

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
  Data <- Data[Data$YEAR == year, ]
  
  DataTA <- DataTA[!is.na(DataTA$AREA), ]
  Data <- Data[!is.na(Data$AREA), ]

  if (nrow(Data) == 0) {
    write(paste("\n----------- check correctness of date by haul - ", year), file = Errors, append = TRUE)
    write("No data found for this year", file = Errors, append = TRUE)
    return(TRUE)
  }

  tf <- as.character(Data$TYPE_OF_FILE[1])
  write(paste("\n----------- check correctness of date by haul in", tf, "-", year), file = Errors, append = TRUE)

  Data$Date <- paste(Data$AREA, Data$COUNTRY, Data$HAUL_NUMBER, Data$DAY, Data$MONTH, Data$YEAR, sep="-")
  TA_df <- DataTA
  TA_df$Date <- paste(TA_df$AREA, TA_df$COUNTRY, TA_df$HAUL_NUMBER, TA_df$DAY, TA_df$MONTH, TA_df$YEAR, sep="-")

  numberError <- 0
  unique_dataset <- Data[!duplicated(Data$Date), ]
  
  for (i in 1:nrow(unique_dataset)) {
    if (!(unique_dataset$Date[i] %in% TA_df$Date)) {
      if (tf != "TL") {
        write(paste("Haul", unique_dataset$HAUL_NUMBER[i], ", code species", unique_dataset$GENUS[i], unique_dataset$SPECIES[i], ": the date is not consistent with the date reported in TA."), file = Errors, append = TRUE)
      } else {
        write(paste("Haul", unique_dataset$HAUL_NUMBER[i], ": the date is not consistent with the date reported in TA."), file = Errors, append = TRUE)
      }
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