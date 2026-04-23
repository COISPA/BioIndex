#' Check numeric range (RoME)
#' @description
#' The function checks whether the values contained in specific fields are consistent with the allowed ranges.
#'
#' @param ResultData data frame in MEDITS tables
#' @param Field field of the table to be checked
#' @param Values vector of the allowed values
#' @param year reference year for the analysis
#' @param wd working directory
#' @param suffix name of the log file
#' @param verbose boolean. If TRUE messages are prompted in the console
#' @export
check_numeric_range <- function(ResultData, Field, Values, year, wd=NA, suffix=NA, verbose=FALSE) {

  if (is.na(wd)) {
    wd <- tempdir()
    if (verbose) message(paste("Missing working directory. Results are saved in the temporary folder: ", wd))
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

  # Filtering data for the selected year
  if (length(year) != 1 || is.na(year)) {
    stop("Argument 'year' must be a single non-NA value")
  }
  df <- ResultData[ResultData$YEAR == year, , drop = FALSE]
  if (nrow(df) == 0) {
    write(paste("\n----------- check range of values for field:", Field, "-", year), file = Errors, append = TRUE)
    write(sprintf("No data found for year %s", year), file = Errors, append = TRUE)
    return(TRUE)
  }
  
  tf <- df$TYPE_OF_FILE[1]
  write(paste("\n----------- check range of values for field:", Field, "-", year), file = Errors, append = TRUE)

  indexcol <- which(names(df) == Field)
  Valuesf <- as.numeric(Values)
  lrange <- length(Valuesf)

  # Remove NAs for range check
  df_sub <- df[!is.na(df[[Field]]), ]
  numberError <- 0

  if (nrow(df_sub) > 0) {
    if (lrange == 2) {
      for (k in 1:nrow(df_sub)) {
        val <- df_sub[k, indexcol]
        if (!(val >= Valuesf[1] & val <= Valuesf[2])) {
          write(paste("Haul", df_sub$HAUL_NUMBER[k], ": value (", val, ") out of allowed range for", Field, "in", tf), file = Errors, append = TRUE)
          numberError <- numberError + 1
        }
      }
    } else if (lrange > 2) {
      for (k in 1:nrow(df_sub)) {
        val <- df_sub[k, indexcol]
        if (!((val >= Valuesf[1] & val <= Valuesf[2]) | val %in% Valuesf[3:lrange])) {
          write(paste("Haul", df_sub$HAUL_NUMBER[k], ": value (", val, ") out of allowed range for", Field, "in", tf), file = Errors, append = TRUE)
          numberError <- numberError + 1
        }
      }
    }
  }

  if (numberError == 0) {
    write(paste("No error occurred for field", Field, "in", tf), file = Errors, append = TRUE)
    return(TRUE)
  } else {
    return(FALSE)
  }
}
