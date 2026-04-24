#' Conversion of decimal degrees coordinates in MEDITS format
#' @description Conversion of decimal degrees coordinates in MEDITS format
#' @param data data frame of the hauls data (TA, table A) in MEDITS format
#' @return The function returns the data frame of the TA (table A) reporting the coordinates in MEDITS format.
#' @export
dd.to.MEDITS <- function(data)  {

  deg <- floor(data$SHOOTING_LATITUDE)
  dec <- data$SHOOTING_LATITUDE - (deg)
  dd  <- (deg*100)+(dec*60)
  data$SHOOTING_LATITUDE <- abs(dd)

  ##
  deg <- floor(data$SHOOTING_LONGITUDE)
  dec <- data$SHOOTING_LONGITUDE - (deg)
  dd  <- (deg*100)+(dec*60)
  data$SHOOTING_LONGITUDE <- abs(dd)

  ##
  deg <- floor(data$HAULING_LATITUDE)
  dec <- data$HAULING_LATITUDE - (deg)
  dd  <- (deg*100)+(dec*60)
  data$HAULING_LATITUDE <- abs(dd)

  ##
  deg <- floor(data$HAULING_LONGITUDE)
  dec <- data$HAULING_LONGITUDE - (deg)
  dd  <- (deg*100)+(dec*60)
  data$HAULING_LONGITUDE <- abs(dd)

  return(data)
}

