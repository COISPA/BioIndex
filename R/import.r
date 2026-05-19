#' @importFrom stats AIC coef coefficients cor.test dchisq deviance
#' @importFrom stats df.residual lm na.exclude nls pnorm predict rnorm sd time quantile
#' @importFrom utils data write.table read.table
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange mutate filter select group_by summarize
#' @importFrom terra unwrap
#' @importFrom graphics box
#' @importFrom grDevices dev.off jpeg heat.colors
#' @importFrom utils data write.table
#' @importFrom shinyjs useShinyjs
NULL

# Solve R CMD check notes on global variables
utils::globalVariables(c("AGE", "LENGTH", "LENGTH_CLASS", "lab_y", "med_bathy", "strata",
                         "N_km2", "kg_km2", "YEAR", "STRATUM_CODE", "SWEPT_AREA",
                         "n_raised", "num_weighted", "LENGTH_CLASSES_CODE",
                         "GENUS", "SPECIES", "COUNTRY", "MEAN_DEPTH",
                         "NB_OF_FEMALES", "NB_OF_MALES", "NB_FM", "Indices_F",
                         "Indices_FM", "sr", "variance", "stratum", "value",
                         "year", "STRATUM", "numb_per_stratum", "id", "ALL.STRATA",
                         "longitude", "latitude", "group", "V1", "V2", "V3"))