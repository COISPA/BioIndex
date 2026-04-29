#' Merge TA–TB and TA–TC tables (MEDITS protocol)
#'
#' @author Walter Zupa \email{zupa@fondazionecoispa.org}
#' @description
#' This function is the primary entry point for integrating MEDITS survey data
#' (e.g. TA, TB, and TC tables), producing unified datasets essential for
#' population analysis and the estimation of demographic indicators.
#' The routine:
#' \itemize{
#'   \item filters hauls by validity and by country;
#'   \item performs a full set of quality checks (times, positions, wing opening,
#'         haul/date consistency) using internal validation functions
#'         (based on \code{RoME} version 0.2.3);
#'   \item merges TA with TB and TA with TC, respectively;
#'   \item computes swept area, mean positions, depth stratum, raising
#'         factors, densities and biomasses;
#'   \item (optionally) writes the two merged tables to
#'         \code{file.path(wd, "output")}.
#' }
#'
#' @param ta      A MEDITS or MEDITS-like \strong{TA} table
#'                (columns listed in \code{BioIndex::TA_cols}).
#' @param tb      A MEDITS or MEDITS-like \strong{TB} table
#'                (columns listed in \code{BioIndex::TB_cols}).
#' @param tc      A MEDITS or MEDITS-like \strong{TC} table
#'                (columns listed in \code{BioIndex::TC_cols}).
#' @param species Character. MEDITS 7-letter RubIN code (e.g. \code{"MERLMER"}):
#'                the first 4 letters are used as \emph{GENUS}, the last 3 as
#'                \emph{SPECIES}.
#' @param country Character vector of MEDITS country codes to keep.
#'                Use \code{"all"} (default) to include every country present
#'                in \code{ta}.
#' @param strata  A data-frame with the depth-stratification scheme adopted
#'                by the MEDITS survey. Defaults to
#'                \code{BioIndex::strata_scheme}.
#' @param wd      Working directory. When \code{save = TRUE}, the merged
#'                tables are written to \code{file.path(wd, "output")}.
#' @param save    Logical. If \code{TRUE} (default) the function writes
#'                \code{"mergeTATB_<GENUS><SPECIES>.csv"} and
#'                \code{"mergeTATC_<GENUS><SPECIES>.csv"}.
#' @param verbose Logical (default \code{TRUE}); prints progress messages.
#'
#' @return A list of two \code{data.frame}s:
#'   \describe{
#'     \item{\code{merge_TA_TB}}{One row per haul with TA metadata, TB catch
#'           totals, depth stratum, densities and biomasses.}
#'     \item{\code{merge_TA_TC}}{One row per haul/length-class/sex/maturity
#'           with TA metadata, TC counts, depth stratum, raising factor,
#'           densities and biomasses.}
#'   }
#'
#' @details
#' The implementation reproduces the official BioIndex workflow.
#'
#' \strong{Data Validation:}
#' The function automatically performs syntactic data validation using internal
#' routines based on \code{RoME} version 0.2.3 logic to ensure data integrity
#' and standardisation.
#'
#' \strong{Optimisations:}
#' Two micro-optimisations are included:
#' \enumerate{
#'   \item vectorised replacement of missing \code{"NA NA"} entries
#'         in both merges;
#'   \item a single loop over depth strata (instead of haul × stratum) to
#'         assign \code{STRATUM_CODE}.
#' }
#' Numerical outputs are identical to the reference routine.
#'
#' @note
#' This version of \code{BioIndex} uses internal validation functions based on
#' \strong{RoME version 0.2.3}. The package is fully self-sufficient and does
#' not require external validation packages to be installed.
#'
#' @examples
#' \donttest{
#' data("TA", package = "BioIndex")
#' data("TB", package = "BioIndex")
#' data("TC", package = "BioIndex")
#' m <- merge_TATBTC(
#'   ta      = TA,
#'   tb      = TB,
#'   tc      = TC,
#'   species = "MERLMER",
#'   country = "all",
#'   wd      = tempdir(),
#'   verbose = TRUE
#' )
#' mTATB <- m[[1]]   # TA–TB merged table
#' mTATC <- m[[2]]   # TA–TC merged table
#' }
#'
#' @export
#' @importFrom hms hms
#' @importFrom utils write.table
merge_TATBTC <- function(ta, tb, tc,
                         species,
                         country = "all",
                         strata  = BioIndex::strata_scheme,
                         wd      = NA,
                         save    = TRUE,
                         verbose = TRUE) {

  if (is.na(wd)) { wd <- tempdir() }

  strata_scheme <- strata  # local alias

  ## ------------------------------------------------------------------ ##
  ##  Column templates                                                  ##
  ## ------------------------------------------------------------------ ##
  TA_cols <- BioIndex::TA_cols
  TB_cols <- BioIndex::TB_cols
  TC_cols <- BioIndex::TC_cols

  ## ------------------------------------------------------------------ ##
  ##  Species code parsing                                              ##
  ## ------------------------------------------------------------------ ##
  species[2] <- substr(species, 5, 7)
  species[1] <- substr(species[1], 1, 4)
  sspp  <- paste(species[1], species[2], sep = "")
  GENERE <- species[1]
  SPECIE <- species[2]

  ## ------------------------------------------------------------------ ##
  ##  Input data frames                                                 ##
  ## ------------------------------------------------------------------ ##
  TA <- ta
  TB <- tb
  TC <- tc

  ## ------------------------------------------------------------------ ##
  ##  COUNTRY filter                                                    ##
  ## ------------------------------------------------------------------ ##
  TA_country <- sort(unique(as.character(TA$COUNTRY)))
  TB_country <- sort(unique(as.character(TB$COUNTRY)))
  TC_country <- sort(unique(as.character(TC$COUNTRY)))
  l_country  <- length(TA_country)

  for (i in seq_len(l_country)) {
    if (!(TA_country[i] %in% TB_country))
      warning(paste("The country ", TA_country[i], " is not present in the TB file.", sep = ""))
    if (!(TA_country[i] %in% TC_country))
      warning(paste("The country ", TA_country[i], " is not present in the TC file.", sep = ""))
  }

  if (l_country == 1) {
    country_analysis <- TA_country
  } else {
    country_analysis <- if (any(country %in% "all")) TA_country else country
  }

  TA <- TA[TA$COUNTRY %in% country_analysis, ]
  TB <- TB[TB$COUNTRY %in% country_analysis, ]
  TC <- TC[TC$COUNTRY %in% country_analysis, ]

  ## ------------------------------------------------------------------ ##
  ##  Unique haul identifiers                                           ##
  ## ------------------------------------------------------------------ ##
  id_TA <- data.frame(id = paste(TA$AREA, TA$COUNTRY, TA$YEAR, "_",
                                 TA$VESSEL, TA$MONTH, TA$DAY, "_",
                                 TA$HAUL_NUMBER, sep = ""))
  id_TB <- data.frame(id = paste(TB$AREA, TB$COUNTRY, TB$YEAR, "_",
                                 TB$VESSEL, TB$MONTH, TB$DAY, "_",
                                 TB$HAUL_NUMBER, sep = ""))
  id_TC <- data.frame(id = paste(TC$AREA, TC$COUNTRY, TC$YEAR, "_",
                                 TC$VESSEL, TC$MONTH, TC$DAY, "_",
                                 TC$HAUL_NUMBER, sep = ""))

  id_invalid <- id_TA$id[which(TA$VALIDITY == "I")]

  ## Data validation: local BioIndex routines (based on RoME 0.2.3)     ##
  suffix <- paste(Sys.Date(), format(Sys.time(), "_time_h%Hm%Ms%OS0"), sep = "")
  for (yy in sort(unique(TA$YEAR))) {
      TAy <- TA[TA$YEAR == yy, ]
      TBy <- TB[TB$YEAR == yy, ]
      TCy <- TC[TC$YEAR == yy, ]

      if (!check_dictionary(TAy, "SHOOTING_TIME", 0:2400, yy, file.path(wd, "output"), suffix))
          stop("SHOOTING_TIME out of range.")
      if (!check_dictionary(TAy, "HAULING_TIME", 0:2400, yy, file.path(wd, "output"), suffix))
          stop("HAULING_TIME out of range.")
      if (!check_dictionary(TAy, "WING_OPENING", c(30, 50:250), yy, file.path(wd, "output"), suffix))
          stop("WING_OPENING out of range.")
      if (!check_numeric_range(TAy, "SHOOTING_LATITUDE", c(3020, 4730), yy, file.path(wd, "output"), suffix))
          stop("SHOOTING_LATITUDE out of range.")
      if (!check_numeric_range(TAy, "HAULING_LATITUDE", c(3020, 4730), yy, file.path(wd, "output"), suffix))
          stop("HAULING_LATITUDE out of range.")
      if (!check_numeric_range(TAy, "SHOOTING_LONGITUDE", c(0, 4200), yy, file.path(wd, "output"), suffix))
          stop("SHOOTING_LONGITUDE out of range.")
      if (!check_numeric_range(TAy, "HAULING_LONGITUDE", c(0, 4200), yy, file.path(wd, "output"), suffix))
          stop("HAULING_LONGITUDE out of range.")
      if (!check_date_haul(TAy, TBy, yy, file.path(wd, "output"), suffix))
          stop("Date in TB not consistent with TA.")
      if (!check_date_haul(TAy, TCy, yy, file.path(wd, "output"), suffix))
          stop("Date in TC not consistent with TA.")
      if (!check_hauls_TBTA(TAy, TBy, yy, file.path(wd, "output"), suffix))
          stop("Hauls in TB not consistent with TA.")
  }

  ## ------------------------------------------------------------------ ##
  ##  Data preparation for merging                                      ##
  ## ------------------------------------------------------------------ ##
  colnames(TA)[colnames(TA) == "AREA"] <- "GSA"
  colnames(TB)[colnames(TB) == "AREA"] <- "GSA"
  colnames(TC)[colnames(TC) == "AREA"] <- "GSA"

  TA_merge <- cbind(id_TA, TA)
  TA_merge <- TA_merge[!TA_merge$id %in% id_invalid, ]
  TA_merge <- TA_merge[, which(colnames(TA_merge) %in% TA_cols)]

  TB_merge <- cbind(id_TB, TB)
  TB_merge$GENUS   <- as.character(TB_merge$GENUS)
  TB_merge$SPECIES <- as.character(TB_merge$SPECIES)
  TB_merge <- TB_merge[!TB_merge$id %in% id_invalid, ]
  TB_merge <- TB_merge[TB_merge$GENUS == species[1] &
                         TB_merge$SPECIES == species[2],
                       which(colnames(TB_merge) %in% TB_cols)]

  TC_merge <- cbind(id_TC, TC)
  TC_merge$GENUS   <- as.character(TC_merge$GENUS)
  TC_merge$SPECIES <- as.character(TC_merge$SPECIES)
  TC_merge <- TC_merge[!TC_merge$id %in% id_invalid, ]
  TC_merge <- TC_merge[TC_merge$GENUS == species[1] &
                         TC_merge$SPECIES == species[2],
                       which(colnames(TC_merge) %in% TC_cols)]

  #######################################################################
  ##                 MERGE TA – TB  (optimised)                        ##
  #######################################################################
  if (verbose) message("- Merging TA-TB files")

  merge_TATB <- merge(TA_merge, TB_merge, by = "id", all.x = TRUE)
  merge_TATB$MEDITS_CODE <- paste(merge_TATB$GENUS, merge_TATB$SPECIES)

  ## -------- Optimization 1: remove row-by-row loop ------------------ ##
  na_idx <- merge_TATB$MEDITS_CODE == "NA NA"
  if (any(na_idx)) {
    merge_TATB$MEDITS_CODE[na_idx]               <- "NA"
    merge_TATB$GENUS[na_idx]                     <- -1
    merge_TATB$SPECIES[na_idx]                   <- -1
    merge_TATB$TOTAL_WEIGHT_IN_THE_HAUL[na_idx]  <- 0
    merge_TATB$TOTAL_NUMBER_IN_THE_HAUL[na_idx]  <- 0
    merge_TATB$NB_OF_FEMALES[na_idx]             <- 0
    merge_TATB$NB_OF_MALES[na_idx]               <- 0
    merge_TATB$NB_OF_UNDETERMINED[na_idx]        <- 0
  }
  ## ------------------------------------------------------------------ ##

  coord <- convert_coordinates(merge_TATB)

  merge_TATB$MEAN_LATITUDE      <- (merge_TATB$SHOOTING_LATITUDE + merge_TATB$HAULING_LATITUDE) / 2
  merge_TATB$MEAN_LATITUDE_DEC  <- (coord$lat_start + coord$lat_end) / 2
  merge_TATB$MEAN_LONGITUDE     <- (merge_TATB$SHOOTING_LONGITUDE + merge_TATB$HAULING_LONGITUDE) / 2
  merge_TATB$MEAN_LONGITUDE_DEC <- (coord$lon_start + coord$lon_end) / 2
  merge_TATB$MEAN_DEPTH         <- (merge_TATB$SHOOTING_DEPTH + merge_TATB$HAULING_DEPTH) / 2
  merge_TATB$SWEPT_AREA         <- merge_TATB$DISTANCE * merge_TATB$WING_OPENING / 10000000

  strata_subset <- strata_scheme[
    strata_scheme$GSA == unique(merge_TATB$GSA) &
      strata_scheme$COUNTRY %in% as.character(unique(merge_TATB$COUNTRY)), ]

  ## -------- Optimization 2: single loop over strata ------------------ ##
  merge_TATB$STRATUM_CODE <- NA
  for (j in seq_len(nrow(strata_subset))) {
    lower <- strata_subset[j, 4]; upper <- strata_subset[j, 5]; code <- strata_subset[j, 3]
    if (code == 1) {
      idx <- floor(merge_TATB$MEAN_DEPTH) >= lower & floor(merge_TATB$MEAN_DEPTH) <= upper
    } else {
      idx <- floor(merge_TATB$MEAN_DEPTH)  > lower & floor(merge_TATB$MEAN_DEPTH) <= upper
    }
    merge_TATB$STRATUM_CODE[idx] <- code
  }
  ## ------------------------------------------------------------------ ##

  ## Remaining TA-TB calculations unchanged ---------------------------- ##
  hour_shooting <- ifelse(nchar(merge_TATB$SHOOTING_TIME) == 4,
                          substr(merge_TATB$SHOOTING_TIME, 1, 2),
                          substr(merge_TATB$SHOOTING_TIME, 1, 1))
  min_shooting  <- ifelse(nchar(merge_TATB$SHOOTING_TIME) == 4,
                          substr(merge_TATB$SHOOTING_TIME, 3, 4),
                          substr(merge_TATB$SHOOTING_TIME, 2, 3))
  hour_hauling  <- ifelse(nchar(merge_TATB$HAULING_TIME) == 4,
                          substr(merge_TATB$HAULING_TIME, 1, 2),
                          substr(merge_TATB$HAULING_TIME, 1, 1))
  min_hauling   <- ifelse(nchar(merge_TATB$HAULING_TIME) == 4,
                          substr(merge_TATB$HAULING_TIME, 3, 4),
                          substr(merge_TATB$HAULING_TIME, 2, 3))

  hms_shooting <- hms(rep(0, length(hour_shooting)),
                      as.numeric(min_shooting), as.numeric(hour_shooting))
  hms_hauling  <- hms(rep(0, length(hour_hauling)),
                      as.numeric(min_hauling),  as.numeric(hour_hauling))
  duration     <- as.numeric(hms_hauling - hms_shooting) / 3600

  merge_TATB$SHOOTING_TIME <- hms_shooting
  merge_TATB$HAULING_TIME  <- hms_hauling
  merge_TATB$N_h    <- merge_TATB$TOTAL_NUMBER_IN_THE_HAUL / duration
  merge_TATB$N_km2  <- merge_TATB$TOTAL_NUMBER_IN_THE_HAUL / merge_TATB$SWEPT_AREA
  merge_TATB$kg_h   <- merge_TATB$TOTAL_WEIGHT_IN_THE_HAUL / duration / 1000
  merge_TATB$kg_km2 <- merge_TATB$TOTAL_WEIGHT_IN_THE_HAUL / merge_TATB$SWEPT_AREA / 1000
  merge_TATB$MEDITS_CODE <- as.character(merge_TATB$MEDITS_CODE)

  if (save) {
    if (!dir.exists(file.path(wd, "output"))) dir.create(file.path(wd, "output"), recursive = TRUE, showWarnings = FALSE)
    write.table(merge_TATB,
                file.path(wd, "output", paste0("mergeTATB_", species[1], species[2], ".csv")),
                sep = ";", row.names = FALSE)
  }
  if (verbose) message("TA-TB files correctly merged")

  #######################################################################
  ##                 MERGE TA – TC  (optimised)                        ##
  #######################################################################
  if (verbose) message("- Merging TA-TC files")

  merge_TATC <- merge(TA_merge, TC_merge, by = "id", all.x = TRUE, all.y = TRUE)
  merge_TATC$MEDITS_CODE <- paste(merge_TATC$GENUS, merge_TATC$SPECIES)

  ## Same “NA NA” optimization ----------------------------------------- ##
  na_idx2 <- merge_TATC$MEDITS_CODE == "NA NA"
  if (any(na_idx2)) {
    merge_TATC$MEDITS_CODE[na_idx2] <- -1
    merge_TATC[na_idx2,
               c("GENUS","SPECIES","SEX","LENGTH_CLASSES_CODE",
                 "WEIGHT_OF_THE_FRACTION","WEIGHT_OF_THE_SAMPLE_MEASURED",
                 "MATURITY","MATSUB","LENGTH_CLASS")] <- -1
    merge_TATC$WEIGHT_OF_THE_FRACTION[na_idx2] <- 0
    merge_TATC$WEIGHT_OF_THE_SAMPLE_MEASURED[na_idx2] <- 0
    merge_TATC$NUMBER_OF_INDIVIDUALS_IN_THE_LENGTH_CLASS_AND_MATURITY_STAGE[na_idx2] <- 0
  }

  coord2 <- convert_coordinates(merge_TATC)

  merge_TATC$MEAN_LATITUDE      <- (merge_TATC$SHOOTING_LATITUDE + merge_TATC$HAULING_LATITUDE) / 2
  merge_TATC$MEAN_LATITUDE_DEC  <- (coord2$lat_start + coord2$lat_end) / 2
  merge_TATC$MEAN_LONGITUDE     <- (merge_TATC$SHOOTING_LONGITUDE + merge_TATC$HAULING_LONGITUDE) / 2
  merge_TATC$MEAN_LONGITUDE_DEC <- (coord2$lon_start + coord2$lon_end) / 2
  merge_TATC$MEAN_DEPTH         <- (merge_TATC$SHOOTING_DEPTH + merge_TATC$HAULING_DEPTH) / 2
  merge_TATC$SWEPT_AREA         <- merge_TATC$DISTANCE * merge_TATC$WING_OPENING / 10000000

  strata_subset2 <- strata_scheme[
    strata_scheme$GSA == unique(merge_TATB$GSA) &
      strata_scheme$COUNTRY %in% as.character(unique(merge_TATB$COUNTRY)), ]

  ## Optimization 2 also on TA-TC -------------------------------------- ##
  merge_TATC$STRATUM_CODE <- NA
  for (j in seq_len(nrow(strata_subset2))) {
    lower <- strata_subset2[j, 4]; upper <- strata_subset2[j, 5]; code <- strata_subset2[j, 3]
    idx <- merge_TATC$MEAN_DEPTH > lower & merge_TATC$MEAN_DEPTH <= upper
    if (code == 1)
      idx <- merge_TATC$MEAN_DEPTH >= lower & merge_TATC$MEAN_DEPTH <= upper
    merge_TATC$STRATUM_CODE[idx] <- code
  }
  ## ------------------------------------------------------------------ ##

  hour_shooting <- ifelse(nchar(merge_TATC$SHOOTING_TIME) == 4,
                          substr(merge_TATC$SHOOTING_TIME, 1, 2),
                          substr(merge_TATC$SHOOTING_TIME, 1, 1))
  min_shooting  <- ifelse(nchar(merge_TATC$SHOOTING_TIME) == 4,
                          substr(merge_TATC$SHOOTING_TIME, 3, 4),
                          substr(merge_TATC$SHOOTING_TIME, 2, 3))
  hour_hauling  <- ifelse(nchar(merge_TATC$HAULING_TIME) == 4,
                          substr(merge_TATC$HAULING_TIME, 1, 2),
                          substr(merge_TATC$HAULING_TIME, 1, 1))
  min_hauling   <- ifelse(nchar(merge_TATC$HAULING_TIME) == 4,
                          substr(merge_TATC$HAULING_TIME, 3, 4),
                          substr(merge_TATC$HAULING_TIME, 2, 3))

  hms_shooting <- hms(rep(0, length(hour_shooting)),
                      as.numeric(min_shooting), as.numeric(hour_shooting))
  hms_hauling  <- hms(rep(0, length(hour_hauling)),
                      as.numeric(min_hauling),  as.numeric(hour_hauling))
  duration     <- as.numeric(hms_hauling - hms_shooting) / 3600

  merge_TATC$SHOOTING_TIME <- hms_shooting
  merge_TATC$HAULING_TIME  <- hms_hauling

  merge_TATC$MATURITY_STAGE <- with(merge_TATC, ifelse(
    MATURITY %in% c(-1, 0, 1), MATURITY,
    ifelse(MATSUB %in% c("A","B","C","D","E","F"),
           paste(MATURITY, MATSUB, sep = ""), MATURITY)))

  merge_TATC$RAISING_FACTOR <- 0
  pos <- merge_TATC$WEIGHT_OF_THE_FRACTION > 0
  merge_TATC$RAISING_FACTOR[pos] <- merge_TATC$WEIGHT_OF_THE_FRACTION[pos] /
    merge_TATC$WEIGHT_OF_THE_SAMPLE_MEASURED[pos]
  merge_TATC$RAISING_FACTOR[merge_TATC$RAISING_FACTOR > 0 &
                              merge_TATC$RAISING_FACTOR < 1] <- 1

  merge_TATC$N_h    <- merge_TATC$NUMBER_OF_INDIVIDUALS_IN_THE_LENGTH_CLASS_AND_MATURITY_STAGE *
    merge_TATC$RAISING_FACTOR / duration
  merge_TATC$N_km2  <- merge_TATC$NUMBER_OF_INDIVIDUALS_IN_THE_LENGTH_CLASS_AND_MATURITY_STAGE *
    merge_TATC$RAISING_FACTOR / merge_TATC$SWEPT_AREA
  merge_TATC$kg_h   <- merge_TATC$WEIGHT_OF_THE_SAMPLE_MEASURED / duration / 1000
  merge_TATC$kg_km2 <- merge_TATC$WEIGHT_OF_THE_SAMPLE_MEASURED / merge_TATC$SWEPT_AREA / 1000

  if (save) {
    if (!dir.exists(file.path(wd, "output"))) dir.create(file.path(wd, "output"), recursive = TRUE, showWarnings = FALSE)
    write.table(merge_TATC,
                file.path(wd, "output", paste0("mergeTATC_", species[1], species[2], ".csv")),
                sep = ";", row.names = FALSE)
  }
  if (verbose) message("TA-TC files correctly merged")

  ## ------------------------------------------------------------------ ##
  ##  Return                                                            ##
  ## ------------------------------------------------------------------ ##
  return(list(merge_TA_TB = merge_TATB,
              merge_TA_TC = merge_TATC))
}
