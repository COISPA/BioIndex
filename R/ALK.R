#' ALK (Age Length Key)
#'
#' @param ta MEDITS or MEDITS-like TA table
#' @param te MEDITS or MEDITS-like TE table
#' @param sp species rubin code (MEDITS format, e.g. "MERLMER")
#' @param GSA reference GSA for the analysis
#' @param country reference country
#' @param nyears number of years of the time series to be considered in the analysis
#' @param wd path of the working directory
#' @param verbose boolean. If TRUE messages are promted in the console
#' @export
#' @import ggplot2
ALK <- function(ta,te,sp,GSA,country="all",nyears=NA,wd=NA,verbose=FALSE){

    if (FALSE) {
        ta
        te

        wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex_appoggio\\test"
        ta <- RoME::TA
        te <- RoME::TE
        sp="MERLMER"
        GSA=10
        country="all"
        nyears=NA
        # wd=NA
        verbose=FALSE

        ALK(ta,te,sp,GSA,country="all",nyears=NA,wd=NA,verbose=FALSE)
    }

    if (is.na(wd)) {
        wd <- getwd()
    }

    TA <- ta
    TE <- te

    if (verbose) {
        cat("\n###########################\n")
        cat("Age-Length keys (ALK)\n")
        cat("###########################\n")
        cat("\n")
    }

      if (nrow(TE)>0 & nrow(TA)>0) {

          if (file.exists(file.path(wd,"output"))) {
              setwd(wd)
          } else {
              dir.create(file.path(wd,"output"), showWarnings = FALSE)
              setwd(wd)
          }

            if (file.exists(file.path(wd,"output/ALK"))) {
                setwd(wd)
            } else {
                dir.create(file.path(wd,"output/ALK"), showWarnings = FALSE)
                setwd(wd)
            }


            GEAR <- unique(TA$GEAR)[1]
            YEARS <- sort(unique(TE$YEAR))
            nYEARS <- length(YEARS)

            if (is.na(nyears)) {
                n_years <- nYEARS
                if (verbose) { message(paste0("Only ",n_years," have been used in the analysis."))}
            } else if (nyears >= nYEARS){
                n_years <- nYEARS
                if (verbose) { message(paste0("The available time series is shorter than ",nYEARS," years. Only ",n_years," have been used in the analysis."))}
            }

            #YEARS
            if (!is.na(n_years) & n_years <= nYEARS & n_years >0) {
                n_years <- n_years
            } else {
                n_years <- nYEARS
            }

            # SPECIES

                species <- substr(sp[1],1,4)
                species[2] <- substr(sp,5,7)
                sspp <- paste(species[1],species[2], sep="")


            # COUNTRY

                #########################
                ###   filtro COUNTRY  ###
                #########################
                TA_country <- sort(unique(as.character(TA$COUNTRY)))
                TE_country <- sort(unique(as.character(TE$COUNTRY)))

                l_country <- length(TA_country)

                i=1
                for (i in 1:l_country){
                    if (!(TA_country[i] %in% TE_country)){warning(paste("The country ",TA_country[i]," is not present in the TE file.", sep=""))}
                }

                if (l_country==1){
                    check_country ="Y"
                    country_analysis <- TA_country
                } else {

                    if (any(country %in% "all")) {
                        country_analysis <<- TA_country
                    } else {
                        country_analysis <<- country
                    }
                } # close -->if  length(l_country)==1

                TA <- TA[TA$COUNTRY %in% country_analysis , ]
                TE <- TE[TE$COUNTRY %in% country_analysis , ]




                TE <- TE[TE$AREA %in% GSA, ]
                if (nrow(TE)>0){
                       ALKf(TE, GEAR, GSA, country = country_analysis, sp=sspp, years = n_years)
                    }
        } # close (nrow(TE)>0 & nrow(TA)>0)

} # close function
