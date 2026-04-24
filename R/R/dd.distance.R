#' Estimate hauls distances (decimal degrees)
#' @description
#' Function to estimate the hauls length using TA (table A, hauls data) with coordinates in the decimal degrees format (dd.ddd). The distances could be returned expressed in meters, kilometers and nautical miles.
#'
#' @param data data frame of the hauls data (TA, table A) with coordinates reported as decimal degrees
#' @param unit string value indicating the measure unit of the distance. Allowed values: "m" for meters, "km" for kilometers and "NM" for nautical miles.
#' @param verbose give verbose output reporting in the output the selected measure unit of the distance.
#' @return The function returns the vector of the distances expressed in the selected measure unit.
#' @export
dd.distance<-function(data, unit = "m", verbose=TRUE)  {
  N1  <- (((data$SHOOTING_LATITUDE/2)+45)*pi )/180
  N2  <- (((data$HAULING_LATITUDE/2)+45)*pi )/180
  N3  <- atan((pi*(data$HAULING_LONGITUDE-data$SHOOTING_LONGITUDE))/(180*(log(tan(N2))-log(tan(N1)))))
  dist <- abs(60*(data$HAULING_LATITUDE-data$SHOOTING_LATITUDE)/cos(N3))*1852
  if (unit =="m"){
    dist = dist
    if (verbose){
    message("meters")
      }
  }
  if (unit =="km"){
    dist = dist/1000
    if (verbose){
    message("kilometers")
    }
    }
  if (unit =="NM"){
    dist = dist/1852
    if (verbose){
      message("Nautical Miles")
    }
      }
  return(dist)
}
