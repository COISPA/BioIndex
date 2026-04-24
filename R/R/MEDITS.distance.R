#' Estimation of haul distance
#'
#' @param data data frame containing the hauls data (TA, table A).
#' @param unit string value indicating the measure unit of the distance. Allowed values: "m" for meters, "km" for kilometers and "NM" for nautical miles.
#' @param verbose give verbose output reporting in the output the selected measure unit of the distance.
#' @return The function returns the vector of the distances expressed in the selected measure unit.
#' @export
MEDITS.distance<-function(data, unit = "m", verbose = TRUE)  {
  lat_start=data$SHOOTING_LATITUDE
  lon_start= data$SHOOTING_LONGITUDE
  lat_end=data$HAULING_LATITUDE
  lon_end= data$HAULING_LONGITUDE
  LatStartDeg = floor(floor(lat_start)/100);
  LonStartDeg = floor(floor(lon_start)/100);
  LatStartMin=(lat_start-LatStartDeg*100)/60
  LonStartMin=(lon_start-LonStartDeg*100)/60
  LatEndDeg = floor(floor(lat_end)/100);
  LonEndDeg = floor(floor(lon_end)/100);
  LatEndMin=(lat_end-LatEndDeg*100)/60
  LonEndMin=(lon_end-LonEndDeg*100)/60

  lat_start2= LatStartDeg + LatStartMin
  lon_start2 = LonStartDeg + LonStartMin
  lat_end2 = LatEndDeg + LatEndMin
  lon_end2 = LonEndDeg + LonEndMin
  data$SHOOTING_LATITUDE = lat_start2
  data$SHOOTING_LONGITUDE = lon_start2
  data$HAULING_LATITUDE = lat_end2
  data$HAULING_LONGITUDE = lon_end2

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

