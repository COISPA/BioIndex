#' Merge TA and TB tables (haul-level catches)
#'
#' @description
#' Combines haul-level operational data (TA) with species catch data (TB),
#' providing the essential foundation for calculating stratified abundance (N)
#' and biomass (kg) indices.
#' The function keeps only valid hauls, performs quality checks (times,
#' positions, wing opening, haul/date consistency) using internal validation
#' functions (based on \code{RoME} version 0.2.3), computes spatial means,
#' swept area, depth stratum, biomass/density indicators, and—optionally—writes
#' the merged data to disk.
#'
#' @param ta   A `data.frame`/`data.table` containing a full TA dataset
#'             (columns listed in `BioIndex::TA_cols`).
#' @param tb   A `data.frame`/`data.table` containing a full TB dataset
#'             (columns listed in `BioIndex::TB_cols`).
#' @param species Character, MEDITS 7-letter code (e.g. `"ARISFOL"`).
#'                The first 4 characters are interpreted as *GENUS* and
#'                the last 3 as *SPECIES*.
#' @param country Character vector of ISO-3 country codes to keep; use
#'                `"all"` (default) to keep every country present in `ta`.
#' @param strata  Depth-stratum reference table; default is
#'                `BioIndex::strata_scheme`.
#' @param wd      Working directory.  If `NA` (the default when
#'                `save     = FALSE`) no files are written.
#' @param save    Logical.  If `TRUE`, writes
#'                `"mergeTATB_<GENUS><SPECIES>.csv"` to
#'                `file.path(wd, "output")`.
#' @param verbose Logical; if `TRUE` (default) prints progress messages.
#'
#' @return A `data.frame` with one row per haul and the following groups
#'         of variables:
#'         * **TA** metadata (haul position, times, depths, etc.)
#'         * **TB** catch totals (numbers and weight)
#'         * Calculated fields: mean lat/lon, swept area, depth stratum,
#'           density (`N_h`, `N_km2`) and biomass (`kg_h`, `kg_km2`)
#' @details
#' The implementation mirrors the original BioIndex routine.
#'
#' \strong{Data Validation:}
#' The function performs syntactic data validation using internal
#' implementations of the validation routines (based on \code{RoME} v0.2.3).
#'
#' \strong{Optimisations:}
#' Includes two speed-ups: (1) vectorised replacement of missing `"NA NA"` records,
#' (2) a single loop over depth strata instead of a nested haul × stratum
#' loop.  Results are identical to the reference routine.
#'
#' @note
#' This version of \code{BioIndex} uses internal validation functions based on
#' \strong{RoME version 0.2.3}. The package is fully self-sufficient and does
#' not require external validation packages to be installed.
#'
#' @seealso [BioIndex] package documentation.
#' @importFrom hms hms
#' @importFrom utils write.table
#' @examples
#' \donttest{
#' data("TA", package = "BioIndex")
#' data("TB", package = "BioIndex")
#' m_tb <- merge_TATB(
#'   ta       = TA,
#'   tb       = TB,
#'   species  = "ARISFOL",
#'   country  = "ITA",
#'   wd       = tempdir(),
#'   save     = FALSE
#' )
#' head(m_tb)
#' }
#' @export

merge_TATB <- function(ta, tb,
                       species,
                       country = "all",
                       strata  = BioIndex::strata_scheme,
                       wd      = NA,
                       save     = FALSE,
                       verbose = TRUE) {

    ## Working directory check (unchanged) ------------------------------- ##
    if (is.na(wd) & save) {
        save <- FALSE
        if (verbose)
            message("Missing working directory. Results will not be saved locally.")
    }

    strata_scheme <- strata           # local alias

    ## Column templates (unchanged) -------------------------------------- ##
    TA_cols <- BioIndex::TA_cols
    TB_cols <- BioIndex::TB_cols

    ## Species code parsing (unchanged) ---------------------------------- ##
    species[2] <- substr(species, 5, 7)
    species[1] <- substr(species[1], 1, 4)

    ## Input data (unchanged) -------------------------------------------- ##
    TA <- ta
    TB <- tb

    ## COUNTRY filter (unchanged, minus TC notes) ------------------------ ##
    TA_country <- sort(unique(TA$COUNTRY))
    TB_country <- sort(unique(TB$COUNTRY))
    for (cn in TA_country)
        if (!(cn %in% TB_country))
            warning("Country ", cn, " not present in TB file.")

    country_analysis <- if (length(TA_country) == 1) TA_country else
        if (any(country %in% "all")) TA_country else country
    TA <- TA[TA$COUNTRY %in% country_analysis, ]
    TB <- TB[TB$COUNTRY %in% country_analysis, ]

    ## Unique haul IDs (unchanged) --------------------------------------- ##
    uid <- function(df)
        paste(df$AREA, df$COUNTRY, df$YEAR, "_",
              df$VESSEL, df$MONTH, df$DAY, "_",
              df$HAUL_NUMBER, sep = "")
    TA$id <- uid(TA); TB$id <- uid(TB)
    invalid_ids <- TA$id[TA$VALIDITY == "I"]

    ## Data validation: local BioIndex routines (based on RoME 0.2.3)     ##
    suffix <- paste(Sys.Date(), format(Sys.time(), "_time_h%Hm%Ms%OS0"), sep = "")
    for (yy in sort(unique(TA$YEAR))) {
        TAy <- TA[TA$YEAR == yy, ]; TBy <- TB[TB$YEAR == yy, ]
        if (!check_dictionary(TAy,"SHOOTING_TIME",0:2400,yy,paste0(wd,"/output/"),suffix))
            stop("SHOOTING_TIME out of range.")
        if (!check_dictionary(TAy,"HAULING_TIME",0:2400,yy,paste0(wd,"/output/"),suffix))
            stop("HAULING_TIME out of range.")
        if (!check_dictionary(TAy,"WING_OPENING",c(30,50:250),yy,paste0(wd,"/output/"),suffix))
            stop("WING_OPENING out of range.")
        if (!check_numeric_range(TAy,"SHOOTING_LATITUDE",c(3020,4730),yy,paste0(wd,"/output/"),suffix))
            stop("SHOOTING_LATITUDE out of range.")
        if (!check_numeric_range(TAy,"HAULING_LATITUDE",c(3020,4730),yy,paste0(wd,"/output/"),suffix))
            stop("HAULING_LATITUDE out of range.")
        if (!check_numeric_range(TAy,"SHOOTING_LONGITUDE",c(0,4200),yy,paste0(wd,"/output/"),suffix))
            stop("SHOOTING_LONGITUDE out of range.")
        if (!check_numeric_range(TAy,"HAULING_LONGITUDE",c(0,4200),yy,paste0(wd,"/output/"),suffix))
            stop("HAULING_LONGITUDE out of range.")
        if (!check_date_haul(TAy,TBy,yy,paste0(wd,"/output/"),suffix))
            stop("Date in TB not consistent with TA.")
        if (!check_hauls_TBTA(TAy,TBy,yy,paste0(wd,"/output/"),suffix))
            stop("Hauls in TB not consistent with TA.")
    }

    ## Preparation for merge (unchanged) --------------------------------- ##
    names(TA)[names(TA)=="AREA"] <- "GSA"
    names(TB)[names(TB)=="AREA"] <- "GSA"

    TA_ <- TA[!TA$id %in% invalid_ids, TA_cols]
    TB_ <- TB[!TB$id %in% invalid_ids &
                  TB$GENUS == species[1] & TB$SPECIES == species[2], TB_cols]

    ## MERGE TA – TB (unchanged logic, incl. optimisations) -------------- ##
    if (verbose) message("- Merging TA-TB files")
    merge_TATB <- merge(TA_, TB_, by = "id", all.x = TRUE)
    merge_TATB$MEDITS_CODE <- paste(merge_TATB$GENUS, merge_TATB$SPECIES)

    bad <- merge_TATB$MEDITS_CODE == "NA NA"
    if (any(bad)) {
        merge_TATB$MEDITS_CODE[bad]            <- "NA"
        merge_TATB$GENUS[bad]                  <- -1
        merge_TATB$SPECIES[bad]                <- -1
        merge_TATB$TOTAL_WEIGHT_IN_THE_HAUL[bad]  <- 0
        merge_TATB$TOTAL_NUMBER_IN_THE_HAUL[bad]  <- 0
        merge_TATB$NB_OF_FEMALES[bad]             <- 0
        merge_TATB$NB_OF_MALES[bad]               <- 0
        merge_TATB$NB_OF_UNDETERMINED[bad]        <- 0
    }

    coord <- convert_coordinates(merge_TATB)
    merge_TATB$MEAN_DEPTH <- (merge_TATB$SHOOTING_DEPTH + merge_TATB$HAULING_DEPTH) / 2
    merge_TATB$SWEPT_AREA <- merge_TATB$DISTANCE * merge_TATB$WING_OPENING / 1e7

    strata_sub <- strata_scheme[
        strata_scheme$GSA == unique(merge_TATB$GSA) &
            strata_scheme$COUNTRY %in% unique(merge_TATB$COUNTRY), ]
    merge_TATB$STRATUM_CODE <- NA
    for (k in seq_len(nrow(strata_sub))) {
        lo <- strata_sub[k,4]; hi <- strata_sub[k,5]; cd <- strata_sub[k,3]
        sel <- if (cd==1)
            floor(merge_TATB$MEAN_DEPTH) >= lo & floor(merge_TATB$MEAN_DEPTH) <= hi
        else
            floor(merge_TATB$MEAN_DEPTH)  > lo & floor(merge_TATB$MEAN_DEPTH) <= hi
        merge_TATB$STRATUM_CODE[sel] <- cd
    }

    ## Duration & density (original lines preserved) --------------------- ##
    hh  <- ifelse(nchar(merge_TATB$SHOOTING_TIME)==4,
                  substr(merge_TATB$SHOOTING_TIME,1,2),
                  substr(merge_TATB$SHOOTING_TIME,1,1))
    mm  <- ifelse(nchar(merge_TATB$SHOOTING_TIME)==4,
                  substr(merge_TATB$SHOOTING_TIME,3,4),
                  substr(merge_TATB$SHOOTING_TIME,2,3))
    hh2 <- ifelse(nchar(merge_TATB$HAULING_TIME)==4,
                  substr(merge_TATB$HAULING_TIME,1,2),
                  substr(merge_TATB$HAULING_TIME,1,1))
    mm2 <- ifelse(nchar(merge_TATB$HAULING_TIME)==4,
                  substr(merge_TATB$HAULING_TIME,3,4),
                  substr(merge_TATB$HAULING_TIME,2,3))

    h1 <- hms(rep(0, length(hh)),  as.numeric(mm),  as.numeric(hh))
    h2 <- hms(rep(0, length(hh2)), as.numeric(mm2), as.numeric(hh2))
    dur <- as.numeric(h2 - h1) / 3600

    merge_TATB$N_h    <- merge_TATB$TOTAL_NUMBER_IN_THE_HAUL / dur
    merge_TATB$N_km2  <- merge_TATB$TOTAL_NUMBER_IN_THE_HAUL / merge_TATB$SWEPT_AREA
    merge_TATB$kg_h   <- merge_TATB$TOTAL_WEIGHT_IN_THE_HAUL / 1000 / dur
    merge_TATB$kg_km2 <- merge_TATB$TOTAL_WEIGHT_IN_THE_HAUL / 1000 / merge_TATB$SWEPT_AREA

    ## Save & return (unchanged) ----------------------------------------- ##
    if (save)
        write.table(merge_TATB,
                    file = file.path(wd, "output",
                                     paste0("mergeTATB_", species[1], species[2], ".csv")),
                    sep = ";", row.names = FALSE)
    if (verbose) message("TA-TB merge completed")

    merge_TATB
}
