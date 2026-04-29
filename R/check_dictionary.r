#' Check dictionary (RoME)
#' @description
#' The function checks whether the values contained in specific fields are consistent with the allowed values of the dictionaries.
#'
#' @param ResultData data frame in MEDITS tables
#' @param Field field of the table to be checked
#' @param Values vector of the allowed values
#' @param year reference year for the analysis
#' @param wd working directory
#' @param suffix name of the log file
#' @param verbose boolean. If TRUE messages are prompted in the console
#' @note
#' This function is an internal routine based on \strong{RoME version 0.2.3}.
#' It is provided within \code{BioIndex} to ensure the package remains
#' functional and self-sufficient for data validation.
#' @return \code{TRUE} if the validation passes, \code{FALSE} otherwise.
#' @export
check_dictionary <- function(ResultData, Field, Values, year, wd=NA, suffix=NA, verbose=FALSE) {

  if (is.na(wd)) {
    wd <- tempdir()
    if (verbose) message(paste("Missing working directory. Results are saved in the temporary folder: ", wd))
  }

  log_dir <- file.path(wd, "Logfiles")
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  
  if (missing(suffix) || is.na(suffix)) {
    suffix <- paste0(Sys.Date(), format(Sys.time(), "_time_h%Hm%Ms%OS0"))
  }
  
  log_file <- file.path(log_dir, paste0("Logfile_", suffix, ".dat"))
  if (!file.exists(log_file)) file.create(log_file)

  write(paste0("\n----------- check dictionary for field: ", Field, " - ", year), file = log_file, append = TRUE)

  if (length(year) != 1 || is.na(year)) {
    stop("Argument 'year' must be a single non-NA value")
  }
  
  df <- ResultData[ResultData$YEAR == year, , drop = FALSE]
  if (nrow(df) == 0) {
     write(sprintf("No data found for year %s", year), file = log_file, append = TRUE)
     return(TRUE)
  }
  tf <- df$TYPE_OF_FILE[1]

  all_vals <- as.character(df[[Field]])
  Values_char <- as.character(Values)
  allowsNA <- any(is.na(Values))
  allowed_vals <- setdiff(Values_char, NA_character_)

  msgs_na <- character()
  if (!allowsNA) {
    na_idx <- which(is.na(all_vals))
    if (length(na_idx) > 0) {
      msgs_na <- sprintf("Haul %s: value not allowed for %s in %s (actual: NA)", df$HAUL_NUMBER[na_idx], Field, tf)
    }
  }

  empty_idx <- which(!is.na(all_vals) & all_vals == "")
  msgs_empty <- character()
  if (length(empty_idx) > 0) {
    msgs_empty <- switch(tf,
                         TA = sprintf("Haul %s: the field %s is empty in %s", df$HAUL_NUMBER[empty_idx], Field, tf),
                         TB = sprintf("Haul %s %s %s: the field %s is empty in %s", df$HAUL_NUMBER[empty_idx], df$GENUS[empty_idx], df$SPECIES[empty_idx], Field, tf),
                         TC = sprintf("Haul %s %s %s %s %s: the field %s is empty in %s", df$HAUL_NUMBER[empty_idx], df$GENUS[empty_idx], df$SPECIES[empty_idx], df$SEX[empty_idx], df$LENGTH_CLASS[empty_idx], Field, tf),
                         TE = sprintf("Haul %s %s %s %s %s: the field %s is empty in %s", df$HAUL_NUMBER[empty_idx], df$GENUS[empty_idx], df$SPECIES[empty_idx], df$SEX[empty_idx], df$LENGTH_CLASS[empty_idx], Field, tf),
                         TL = sprintf("Haul %s %s %s %s %s: the field %s is empty in %s", df$HAUL_NUMBER[empty_idx], df$GENUS[empty_idx], df$SPECIES[empty_idx], df$SEX[empty_idx], df$LENGTH_CLASS[empty_idx], Field, tf),
                         character()
    )
  }

  inv_idx <- which(!is.na(all_vals) & all_vals != "" & !(all_vals %in% allowed_vals))
  if (Field == "SEX") inv_idx <- inv_idx[all_vals[inv_idx] != "FALSE"]
  msgs_inv <- character()
  if (length(inv_idx) > 0) {
    msgs_inv <- switch(tf,
                       TA = sprintf("Haul %s: value '%s' not allowed for %s in %s", df$HAUL_NUMBER[inv_idx], all_vals[inv_idx], Field, tf),
                       TB = sprintf("Haul %s %s %s: value '%s' not allowed for %s in %s", df$HAUL_NUMBER[inv_idx], df$GENUS[inv_idx], df$SPECIES[inv_idx], all_vals[inv_idx], Field, tf),
                       TC = sprintf("Haul %s %s %s %s %s: value '%s' not allowed for %s in %s", df$HAUL_NUMBER[inv_idx], df$GENUS[inv_idx], df$SPECIES[inv_idx], df$SEX[inv_idx], df$LENGTH_CLASS[inv_idx], all_vals[inv_idx], Field, tf),
                       TE = sprintf("Haul %s %s %s %s %s: value '%s' not allowed for %s in %s", df$HAUL_NUMBER[inv_idx], df$GENUS[inv_idx], df$SPECIES[inv_idx], df$SEX[inv_idx], df$LENGTH_CLASS[inv_idx], all_vals[inv_idx], Field, tf),
                       TL = sprintf("Haul %s: value '%s' not allowed for %s in %s", df$HAUL_NUMBER[inv_idx], all_vals[inv_idx], Field, tf),
                       character()
    )
  }

  msgs <- c(msgs_na, msgs_empty, msgs_inv)
  if (length(msgs) == 0) {
    write(sprintf("No error occurred for field %s in %s", Field, tf), file = log_file, append = TRUE)
    return(TRUE)
  } else {
    write(msgs, file = log_file, append = TRUE)
    return(FALSE)
  }
}
