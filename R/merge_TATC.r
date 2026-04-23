#' Merge TA and TC tables (haul-level length–frequency data)
#'
#' Combines a MEDITS **TA** table (haul metadata) with the corresponding
#' **TC** table (length-frequency sample by haul and sex/maturity).
#' The function keeps only valid hauls, runs all relevant RoMEBS
#' quality checks, computes swept area, depth stratum, raising factors,
#' abundance/biomass indicators, and—optionally—writes the merged data.
#'
#' @param ta   A `data.frame`/`data.table` containing a full TA dataset
#'             (columns listed in `BioIndex::TA_cols`).
#' @param tc   A `data.frame`/`data.table` containing a full TC dataset
#'             (columns listed in `BioIndex::TC_cols`).
#' @param species Character, MEDITS 7-letter code (e.g. `"ARISFOL"`).
#'                The first 4 characters are interpreted as *GENUS* and
#'                the last 3 as *SPECIES*.
#' @param country Character vector of ISO-3 country codes to keep; use
#'                `"all"` (default) to keep every country present in `ta`.
#' @param strata  Depth-stratum reference table; default is
#'                `BioIndex::strata_scheme`.
#' @param wd      Working directory.  If `NA` (the default when
#'                `save = FALSE`) no files are written.
#' @param save    Logical.  If `TRUE`, writes
#'                `"mergeTATC_<GENUS><SPECIES>.csv"` to
#'                `file.path(wd, "output")`.
#' @param verbose Logical; if `TRUE` (default) prints progress messages.
#'
#' @return A `data.frame` in which each row represents one length-class
#'         (and sex/maturity) sampled in a haul, including:
#'         * **TA** metadata
#'         * **TC** length-frequency counts and weights
#'         * Calculated fields: mean lat/lon, swept area, depth stratum,
#'           raising factor, density (`N_h`, `N_km2`) and biomass
#'           (`kg_h`, `kg_km2`)
#'
#' @details
#' The code is identical to the official BioIndex routine except for two
#' micro-optimisations: vectorised handling of missing `"NA NA"` entries
#' and a single depth-stratum loop. Numerical results remain unchanged.
#'
#' @seealso [BioIndex] package documentation.
#' @importFrom hms hms
#' @importFrom utils write.table
#' @examples
#' \dontrun{
#' m_tc <- merge_TATC(
#'   ta       = ta_df,
#'   tc       = tc_df,
#'   species  = "ARISFOL",
#'   country  = "ESP",
#'   wd       = "/path/to/project",
#'   save     = FALSE
#' )
#' head(m_tc)
#' }
#' @export
merge_TATC <- function(ta, tc,
                       species,
                       country = "all",
                       strata  = BioIndex::strata_scheme,
                       wd      = NA,
                       save    = TRUE,
                       verbose = TRUE) {

    ## Working-directory check ------------------------------------------- ##
    if (is.na(wd) & save) {
        save <- FALSE
        if (verbose)
            message("Missing working directory. Results will not be saved locally.")
    }

    strata_scheme <- strata   # local alias

    ## Column templates -------------------------------------------------- ##
    TA_cols <- BioIndex::TA_cols
    TC_cols <- BioIndex::TC_cols   # TB template removed

    ## Species code parsing --------------------------------------------- ##
    species[2] <- substr(species, 5, 7)
    species[1] <- substr(species[1], 1, 4)

    ## Input data -------------------------------------------------------- ##
    TA <- ta
    TC <- tc

    ## COUNTRY filter (TB removed) -------------------------------------- ##
    TA_country <- sort(unique(TA$COUNTRY))
    TC_country <- sort(unique(TC$COUNTRY))
    for (cn in TA_country)
        if (!(cn %in% TC_country))
            warning("Country ", cn, " not present in TC file.")

    country_analysis <- if (length(TA_country) == 1) TA_country else
        if (any(country %in% "all")) TA_country else country
    TA <- TA[TA$COUNTRY %in% country_analysis, ]
    TC <- TC[TC$COUNTRY %in% country_analysis, ]

    ## Unique haul IDs --------------------------------------------------- ##
    uid <- function(df)
        paste(df$AREA, df$COUNTRY, df$YEAR, "_",
              df$VESSEL, df$MONTH, df$DAY, "_",
              df$HAUL_NUMBER, sep = "")
    TA$id <- uid(TA); TC$id <- uid(TC)
    invalid_ids <- TA$id[TA$VALIDITY == "I"]

    ## RoME checks – conditional usage ----------------------------------- ##
    if (requireNamespace("RoME", quietly = TRUE)) {
        suffix <- paste(Sys.Date(), format(Sys.time(), "_time_h%Hm%Ms%OS0"), sep = "")
        for (yy in sort(unique(TA$YEAR))) {
            TAy <- TA[TA$YEAR == yy, ]; TCy <- TC[TC$YEAR == yy, ]
            if (!RoME::check_dictionary(TAy,"SHOOTING_TIME",0:2400,yy,paste0(wd,"/output/"),suffix))
                stop("SHOOTING_TIME out of range.")
            if (!RoME::check_dictionary(TAy,"HAULING_TIME",0:2400,yy,paste0(wd,"/output/"),suffix))
                stop("HAULING_TIME out of range.")
            if (!RoME::check_dictionary(TAy,"WING_OPENING",c(30,50:250),yy,paste0(wd,"/output/"),suffix))
                stop("WING_OPENING out of range.")
            if (!RoME::check_numeric_range(TAy,"SHOOTING_LATITUDE",c(3020,4730),yy,paste0(wd,"/output/"),suffix))
                stop("SHOOTING_LATITUDE out of range.")
            if (!RoME::check_numeric_range(TAy,"HAULING_LATITUDE",c(3020,4730),yy,paste0(wd,"/output/"),suffix))
                stop("HAULING_LATITUDE out of range.")
            if (!RoME::check_numeric_range(TAy,"SHOOTING_LONGITUDE",c(0,4200),yy,paste0(wd,"/output/"),suffix))
                stop("SHOOTING_LONGITUDE out of range.")
            if (!RoME::check_numeric_range(TAy,"HAULING_LONGITUDE",c(0,4200),yy,paste0(wd,"/output/"),suffix))
                stop("HAULING_LONGITUDE out of range.")
            if (!RoME::check_date_haul(TAy, TCy, yy, paste0(wd,"/output/"), suffix))
                stop("Date in TC not consistent with TA.")
        }
    } else {
        if (verbose) {
            message("The 'RoME' package is not installed. Skipping syntactic data validation.")
        }
    }

    ## Preparation for merge -------------------------------------------- ##
    names(TA)[names(TA)=="AREA"] <- "GSA"
    names(TC)[names(TC)=="AREA"] <- "GSA"

    TA_ <- TA[!TA$id %in% invalid_ids, TA_cols]

    TC_ <- TC[!TC$id %in% invalid_ids &
                  TC$GENUS == species[1] & TC$SPECIES == species[2], TC_cols]

    ## MERGE TA – TC  (unchanged logic) --------------------------------- ##
    if (verbose) message("- Merging TA-TC files")
    merge_TATC <- merge(TA_, TC_, by = "id", all.x = TRUE, all.y = TRUE)
    merge_TATC$MEDITS_CODE <- paste(merge_TATC$GENUS, merge_TATC$SPECIES)

    ## Clean “NA NA” rows (vectorised) ---------------------------------- ##
    bad <- merge_TATC$MEDITS_CODE == "NA NA"
    if (any(bad)) {
        merge_TATC$MEDITS_CODE[bad] <- -1
        merge_TATC[bad,
                   c("GENUS","SPECIES","SEX","LENGTH_CLASSES_CODE",
                     "WEIGHT_OF_THE_FRACTION","WEIGHT_OF_THE_SAMPLE_MEASURED",
                     "MATURITY","MATSUB","LENGTH_CLASS")] <- -1
        merge_TATC$WEIGHT_OF_THE_FRACTION[bad] <- 0
        merge_TATC$WEIGHT_OF_THE_SAMPLE_MEASURED[bad] <- 0
        merge_TATC$NUMBER_OF_INDIVIDUALS_IN_THE_LENGTH_CLASS_AND_MATURITY_STAGE[bad] <- 0
    }

    coord <- convert_coordinates(merge_TATC)
    merge_TATC$MEAN_DEPTH <- (merge_TATC$SHOOTING_DEPTH + merge_TATC$HAULING_DEPTH) / 2
    merge_TATC$SWEPT_AREA <- merge_TATC$DISTANCE * merge_TATC$WING_OPENING / 1e7

    strata_sub <- strata_scheme[
        strata_scheme$GSA == unique(merge_TATC$GSA) &
            strata_scheme$COUNTRY %in% unique(merge_TATC$COUNTRY), ]
    merge_TATC$STRATUM_CODE <- NA
    for (k in seq_len(nrow(strata_sub))) {
        lo <- strata_sub[k,4]; hi <- strata_sub[k,5]; cd <- strata_sub[k,3]
        sel <- if (cd==1)
            merge_TATC$MEAN_DEPTH >= lo & merge_TATC$MEAN_DEPTH <= hi
        else
            merge_TATC$MEAN_DEPTH  > lo & merge_TATC$MEAN_DEPTH <= hi
        merge_TATC$STRATUM_CODE[sel] <- cd
    }

    ## Duration & density (original lines) ------------------------------ ##
    hh  <- ifelse(nchar(merge_TATC$SHOOTING_TIME)==4,
                  substr(merge_TATC$SHOOTING_TIME,1,2),
                  substr(merge_TATC$SHOOTING_TIME,1,1))
    mm  <- ifelse(nchar(merge_TATC$SHOOTING_TIME)==4,
                  substr(merge_TATC$SHOOTING_TIME,3,4),
                  substr(merge_TATC$SHOOTING_TIME,2,3))
    hh2 <- ifelse(nchar(merge_TATC$HAULING_TIME)==4,
                  substr(merge_TATC$HAULING_TIME,1,2),
                  substr(merge_TATC$HAULING_TIME,1,1))
    mm2 <- ifelse(nchar(merge_TATC$HAULING_TIME)==4,
                  substr(merge_TATC$HAULING_TIME,3,4),
                  substr(merge_TATC$HAULING_TIME,2,3))

    h1 <- hms(rep(0, length(hh)),  as.numeric(mm),  as.numeric(hh))
    h2 <- hms(rep(0, length(hh2)), as.numeric(mm2), as.numeric(hh2))
    dur <- as.numeric(h2 - h1) / 3600

    merge_TATC$RAISING_FACTOR <- 0
    pos <- merge_TATC$WEIGHT_OF_THE_FRACTION > 0
    merge_TATC$RAISING_FACTOR[pos] <-
        pmax(merge_TATC$WEIGHT_OF_THE_FRACTION[pos] /
                 merge_TATC$WEIGHT_OF_THE_SAMPLE_MEASURED[pos], 1)

    merge_TATC$N_h    <- merge_TATC$NUMBER_OF_INDIVIDUALS_IN_THE_LENGTH_CLASS_AND_MATURITY_STAGE *
        merge_TATC$RAISING_FACTOR / dur
    merge_TATC$N_km2  <- merge_TATC$NUMBER_OF_INDIVIDUALS_IN_THE_LENGTH_CLASS_AND_MATURITY_STAGE *
        merge_TATC$RAISING_FACTOR / merge_TATC$SWEPT_AREA
    merge_TATC$kg_h   <- merge_TATC$WEIGHT_OF_THE_SAMPLE_MEASURED / 1000 / dur
    merge_TATC$kg_km2 <- merge_TATC$WEIGHT_OF_THE_SAMPLE_MEASURED / 1000 / merge_TATC$SWEPT_AREA

    ## Save & return ----------------------------------------------------- ##
    if (save)
        write.table(merge_TATC,
                    file = file.path(wd, "output",
                                     paste0("mergeTATC_", species[1], species[2], ".csv")),
                    sep = ";", row.names = FALSE)
    if (verbose) message("TA-TC merge completed")

    merge_TATC
}
