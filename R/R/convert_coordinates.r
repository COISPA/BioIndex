#' MEDITS coordinates in decimal degrees
#' @description
#' The function returns the data frame of the TA table with the coordinates expressed as decimal degrees.
#'
#' @param Data data frame of TA table
#'
#' @return the function return the same data frame with the coordinates converted in the decimal degrees format
#' @export convert_coordinates

convert_coordinates<-function(Data)  {

  lat_start=Data$SHOOTING_LATITUDE
  lon_start= Data$SHOOTING_LONGITUDE
  lat_end=Data$HAULING_LATITUDE
  lon_end= Data$HAULING_LONGITUDE
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
  Data$lat_start = lat_start2
  Data$lon_start = lon_start2
  Data$lat_end = lat_end2
  Data$lon_end = lon_end2

  # quadrant 3
  Data[Data$SHOOTING_QUADRANT == 3, "lat_start"] <- -1 * Data[Data$SHOOTING_QUADRANT == 3, "lat_start"]
  Data[Data$HAULING_QUADRANT == 3, "lat_end"] <- -1 * Data[Data$HAULING_QUADRANT == 3, "lat_end"]

  # quadrant 5
  Data[Data$SHOOTING_QUADRANT == 5, "lat_start"] <- -1 * Data[Data$SHOOTING_QUADRANT == 5, "lat_start"]
  Data[Data$HAULING_QUADRANT == 5, "lat_end"] <- -1 * Data[Data$HAULING_QUADRANT == 5, "lat_end"]

  Data[Data$SHOOTING_QUADRANT == 5, "lon_start"] <- -1 * Data[Data$SHOOTING_QUADRANT == 5, "lon_start"]
  Data[Data$HAULING_QUADRANT == 5, "lon_end"] <- -1 * Data[Data$HAULING_QUADRANT == 5, "lon_end"]

  # quadrant 7
  Data[Data$SHOOTING_QUADRANT == 7, "lon_start"] <- -1 * Data[Data$SHOOTING_QUADRANT == 7, "lon_start"]
  Data[Data$HAULING_QUADRANT == 7, "lon_end"] <- -1 * Data[Data$HAULING_QUADRANT == 7, "lon_end"]

  return(Data)
}
