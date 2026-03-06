#' TA table headings
#'
#' @name TA_cols
#' @docType data
#' @author Walter Zupa \email{zupa@fondazionecoispa.it}
#' @keywords TA
"TA_cols"

#' TB table headings
#'
#' @name TB_cols
#' @docType data
#' @author Walter Zupa \email{zupa@fondazionecoispa.it}
#' @keywords TB
"TB_cols"

#' TC table headings
#'
#' @name TC_cols
#' @docType data
#' @author Walter Zupa \email{zupa@fondazionecoispa.it}
#' @keywords TC
"TC_cols"

#' stratification scheme
#'
#' @name strata_scheme
#' @docType data
#' @author Walter Zupa \email{zupa@fondazionecoispa.it}
#' @keywords stratification
#' @description
#' Stratification scheme adopted in the bottom trawl demersal survey (e.g. MEDITS survey).
#'
"strata_scheme"


#' stratification
#'
#' @name stratification
#' @docType data
#' @keywords stratification MEDITS
#' @description
#' Data frame containing the surface area relative to the depth strata adopted in the stratification scheme (strata_scheme) of the demersal surveys (e.g. MEDITS survey).
#'
"stratification"

#' stratification scheme (rapa whelk)
#'
#' @name strata_scheme_rapana
#' @docType data
#' @author Walter Zupa \email{zupa@fondazionecoispa.it}
#' @keywords stratification rapa whelk Black Sea
#' @description
#' Stratification scheme adopted in the rapa whelk survey (e.g. Black Sea beam trawl survey).
#'
"strata_scheme_rapana"


#' stratification (rapa whelk)
#'
#' @name stratification_rapana
#' @docType data
#' @keywords stratification rapa whelk Black Sea
#' @description
#' Data frame containing the surface area relative to the depth strata adopted in the stratification scheme (strata_scheme) of the rapa whelk surveys (e.g. Black Sea beam trawl survey).
#'
"stratification_rapana"

















#' TA table example
#'
#' @name TA
#' @docType data
#' @keywords TA MEDITS
"TA"

#' TB table example
#'
#' @name TB
#' @docType data
#' @keywords TB MEDITS
"TB"

#' TC table example
#'
#' @name TC
#' @docType data
#' @keywords TC MEDITS
"TC"


#' centroidi
#'
#' @name centroidi
#' @docType data
#' @keywords centroidi
"centroidi"

#' cgpmgrid
#'
#' @name cgpmgrid
#' @docType data
#' @keywords cgpmgrid
"cgpmgrid"

#' continent
#'
#' @name continent
#' @docType data
#' @keywords continent
"continent"


#' stratum_0_125
#'
#' @name stratum_0_125
#' @docType data
#' @keywords stratum_0_125
"stratum_0_125"

#' stratum_0_200
#'
#' @name stratum_0_200
#' @docType data
#' @keywords stratum_0_200
"stratum_0_200"


#' stratum_0_35
#'
#' @name stratum_0_35
#' @docType data
#' @keywords stratum_0_35
"stratum_0_35"


#' stratum_0_45
#'
#' @name stratum_0_45
#' @docType data
#' @keywords stratum_0_45
"stratum_0_45"


#' stratum_0_800
#'
#' @name stratum_0_800
#' @docType data
#' @keywords stratum_0_800
"stratum_0_800"


#' stratum_200_800
#'
#' @name stratum_200_800
#' @docType data
#' @keywords stratum_200_800
"stratum_200_800"

#' Mediterranean and Black Sea bathymetry (0–1000 m, bathy object)
#'
#' A precomputed \code{bathy} object containing bathymetric data for the Mediterranean Sea
#' and the Black Sea. Depth values are restricted between 0 and -1000 meters. This dataset
#' was downloaded using \code{\link[marmap]{getNOAA.bathy}} with a resolution of 1 arc-minute,
#' and filtered to remove deeper areas.
#'
#' The object can be used directly with functions from the \pkg{marmap} package, such as
#' \code{plot.bathy()} and \code{get.depth()}.
#'
#' @format An object of class \code{bathy} (a matrix with longitude and latitude as axes, and depth in meters)
#'
#' @details
#' The spatial extent includes:
#' \itemize{
#'   \item Longitude: from -6° to 42°
#'   \item Latitude: from 30° to 47°
#'   \item Depth: from 0 to -1000 meters
#' }
#' Only marine cells within this depth range are retained. Land and deeper areas are set to \code{NA}.
#'
#' @usage data(med_bathy)
#'
#' @source NOAA ETOPO1 via \code{marmap::getNOAA.bathy()}
#' @seealso \code{\link[marmap]{getNOAA.bathy}}, \code{\link[marmap]{plot.bathy}}, \code{\link{bubbleplot_RS_by_hauls}}
#' @keywords datasets bathymetry mediterranean blacksea
"med_bathy"
