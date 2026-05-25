#' ALK (Age Length Key)
#' @description
#' Computes the Age-Length Key (ALK) from biological samples, a crucial
#' parameter for converting length distributions into age classes and
#' analyzing stock growth dynamics.
#'
#' @param ta MEDITS or MEDITS-like TA table
#' @param te MEDITS or MEDITS-like TE table
#' @param sp species RUBIN code (MEDITS format, e.g. "MERLMER")
#' @param GSA reference GSA for the analysis
#' @param country reference country
#' @param nyears number of years of the time series to be considered in the analysis
#' @param wd path of the working directory
#' @param save boolean. If TRUE the outputs are saved in the local folder
#' @param verbose boolean. If TRUE messages are promted in the console
#' @return A \code{data.frame} representing the Age-Length Key.
#' @import ggplot2
ALK <- function(ta, te, sp, GSA, country = "all", nyears = NA, wd = NA, save = TRUE, verbose = FALSE) {
    if (is.na(wd)) {
        wd <- tempdir()
    }

    TA <- ta
    TE <- te

    if (verbose) {
        message("\n###########################")
        message("Age-Length keys (ALK)")
        message("###########################")
        message("")
    }

    if (nrow(TE) > 0 & nrow(TA) > 0) {
        if (!dir.exists(file.path(wd, "output"))) {
            dir.create(file.path(wd, "output"), showWarnings = FALSE)
        }

        if (!dir.exists(file.path(wd, "output", "ALK"))) {
            dir.create(file.path(wd, "output", "ALK"), showWarnings = FALSE, recursive = TRUE)
        }


        GEAR <- unique(TA$GEAR)[1]
        YEARS <- sort(unique(TE$YEAR))
        nYEARS <- length(YEARS)

        if (is.na(nyears)) {
            n_years <- nYEARS
            if (verbose) {
                message(paste0("Only ", n_years, " have been used in the analysis."))
            }
        } else if (nyears >= nYEARS) {
            n_years <- nYEARS
            if (verbose) {
                message(paste0("The available time series is shorter than ", nYEARS, " years. Only ", n_years, " have been used in the analysis."))
            }
        }

        # YEARS
        if (!is.na(n_years) & n_years <= nYEARS & n_years > 0) {
            n_years <- n_years
        } else {
            n_years <- nYEARS
        }

        # SPECIES

        species <- substr(sp[1], 1, 4)
        species[2] <- substr(sp, 5, 7)
        sspp <- paste(species[1], species[2], sep = "")


        # COUNTRY

        #########################
        ###   filtro COUNTRY  ###
        #########################

        TA_country <- as.character(unique(TA$COUNTRY))
        TE_country <- as.character(unique(TE$COUNTRY))
        l_country <- length(TA_country)

        for (i in 1:length(TA_country)) {
            if (!(TA_country[i] %in% TE_country)) {
                warning(paste("The country ", TA_country[i], " is not present in the TE file.", sep = ""))
            }
        }

        if (l_country == 1) {
            check_country <- "Y"
            country_analysis <- TA_country
        } else {
            if (any(country %in% "all")) {
                country_analysis <- TA_country
            } else {
                country_analysis <- country
            }
        } # close -->if  length(l_country)==1

        TA <- TA[TA$COUNTRY %in% country_analysis, ]
        TE <- TE[TE$COUNTRY %in% country_analysis, ]

        TE <- TE[TE$AREA %in% GSA, ]
        if (nrow(TE) > 0) {
            ALKf(te = TE, sp = sspp, GEAR = GEAR, GSA = GSA, country = country_analysis, years = n_years, wd = wd, save = save)
        }
    } # close (nrow(TE)>0 & nrow(TA)>0)
}