#' Aggregate MEDITS data across multiple GSAs
#'
#' Aggregate TA, TB, and TC MEDITS tables across a user-defined set of
#' Geographical Sub-Areas (GSAs), and update the associated stratification
#' objects accordingly.
#'
#' This function filters the input TA, TB, and TC tables to retain only the
#' selected GSAs, concatenates the original \code{AREA} and \code{HAUL_NUMBER}
#' values to create unique haul identifiers across GSAs, and replaces the
#' \code{AREA} field with a single aggregated GSA code obtained by collapsing
#' the selected GSA values into one string. The same aggregation is applied to
#' the \code{stratification} and \code{strata_scheme} tables.
#'
#' @param ta A data frame containing the MEDITS \code{TA} table.
#' @param tb A data frame containing the MEDITS \code{TB} table.
#' @param tc A data frame containing the MEDITS \code{TC} table.
#' @param gsas A numeric or character vector of GSA codes to be aggregated.
#' @param strata_scheme A data frame describing the depth-strata scheme for the
#'   selected GSAs. Defaults to \code{BioIndex::strata_scheme}.
#' @param stratification A data frame containing the stratification information
#'   for the selected GSAs. Defaults to \code{BioIndex::stratification}.
#'
#' @details
#' The function performs three main operations:
#' \enumerate{
#'   \item Checks that all requested GSAs are present in the \code{AREA} column
#'   of the input \code{ta}, \code{tb}, and \code{tc} tables.
#'   \item Filters the three MEDITS tables to the selected GSAs and updates
#'   \code{HAUL_NUMBER} by prefixing it with the original \code{AREA} value,
#'   ensuring uniqueness after aggregation.
#'   \item Filters and updates the \code{stratification} and
#'   \code{strata_scheme} tables so that all selected GSAs are represented under
#'   a single aggregated GSA code.
#' }
#'
#' The aggregated GSA code is created using \code{paste(gsas, collapse = "")}.
#' For example, \code{c(17, 18)} becomes \code{"1718"}.
#'
#' @return
#' A list with five elements:
#' \describe{
#'   \item{[[1]]}{The filtered and aggregated \code{TA} table.}
#'   \item{[[2]]}{The filtered and aggregated \code{TB} table.}
#'   \item{[[3]]}{The filtered and aggregated \code{TC} table.}
#'   \item{[[4]]}{The filtered and updated \code{stratification} table.}
#'   \item{[[5]]}{The filtered and updated \code{strata_scheme} table.}
#' }
#'
#' @examples
#' \dontrun{
#' data("strata_scheme", package = "BioIndex")
#' data("stratification", package = "BioIndex")
#'
#' d <- aggregate_gsas(
#'   ta = ta,
#'   tb = tb,
#'   tc = tc,
#'   gsas = c(17, 18),
#'   strata_scheme = BioIndex::strata_scheme,
#'   stratification = BioIndex::stratification
#' )
#'
#' ta_agg <- d[[1]]
#' tb_agg <- d[[2]]
#' tc_agg <- d[[3]]
#' stratification_agg <- d[[4]]
#' strata_scheme_agg <- d[[5]]
#' }
#'
#' @seealso
#' \code{\link[BioIndex:BioIndex]{strata_scheme}},
#' \code{\link[BioIndex:BioIndex]{stratification}}
#'
#' @export

aggregate_gsas <- function(ta,tb,tc,gsas, strata_scheme = BioIndex::strata_scheme, stratification = BioIndex::stratification) {

    if (FALSE) {
        ta <- read.table("D:/OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L/______ MEDITS DATA __OFFICIAL___/MEDBSsurvey/Demersal/TA_MEDITS_FORMAT_2025.csv", sep=";",header = TRUE)

        tb <- read.table("D:/OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L/______ MEDITS DATA __OFFICIAL___/MEDBSsurvey/Demersal/TB_MEDITS_FORMAT_2025.csv", sep=";",header = TRUE)

        tc <- read.table("D:/OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L/______ MEDITS DATA __OFFICIAL___/MEDBSsurvey/Demersal/TC_MEDITS_FORMAT_2025.csv", sep=";",header = TRUE)

        gsas <- c(17,18)

        strata_scheme = BioIndex::strata_scheme
        stratification = BioIndex::stratification

    }

    try(if(!all(gsas %in% unique(ta$AREA))) stop("\nOne or more GSAs in the selection is not included in the TA table"))

    try(if(!all(gsas %in% unique(tb$AREA))) stop("\nOne or more GSAs in the selection is not included in the TB table"))

    try(if(!all(gsas %in% unique(tc$AREA))) stop("\nOne or more GSAs in the selection is not included in the TC table"))

    ta <- ta[ta$AREA %in% gsas, ]
    tb <- tb[tb$AREA %in% gsas, ]
    tc <- tc[tc$AREA %in% gsas, ]
    ta$HAUL_NUMBER <- paste(ta$AREA,ta$HAUL_NUMBER,sep="_")
    tb$HAUL_NUMBER <- paste(tb$AREA,tb$HAUL_NUMBER,sep="_")
    tc$HAUL_NUMBER <- paste(tc$AREA,tc$HAUL_NUMBER,sep="_")
    ta$AREA <- paste(gsas, collapse="")
    tb$AREA <- paste(gsas, collapse="")
    tc$AREA <- paste(gsas, collapse="")

    stratification <- stratification[stratification$GSA %in% gsas, ]
    strata_scheme <- strata_scheme[strata_scheme$GSA %in% gsas, ]

    stratification$GSA <- paste(gsas, collapse="")
    strata_scheme$GSA <- paste(gsas, collapse="")

    return(list(ta,tb,tc,stratification,strata_scheme))
}
