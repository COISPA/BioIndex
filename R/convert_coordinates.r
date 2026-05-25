#' MEDITS coordinates in decimal degrees
#' @description
#' The function returns the data frame of the TA table with the coordinates expressed as decimal degrees.
#'
#' @param Data data frame of TA table
#'
#' @return the function return the same data frame with the coordinates converted in the decimal degrees format
#' @examples
#' data(TA)
#' convert_coordinates(TA)
#' @export convert_coordinates

convert_coordinates<-function(Data)  {

# [ORIGINALE]:
#   lat_start=Data$SHOOTING_LATITUDE
#   lon_start= Data$SHOOTING_LONGITUDE
#   lat_end=Data$HAULING_LATITUDE
#   lon_end= Data$HAULING_LONGITUDE
#   LatStartDeg = floor(floor(lat_start)/100);
#   LonStartDeg = floor(floor(lon_start)/100);
#   LatStartMin=(lat_start-LatStartDeg*100)/60
#   LonStartMin=(lon_start-LonStartDeg*100)/60
#   LatEndDeg = floor(floor(lat_end)/100);
#   LonEndDeg = floor(floor(lon_end)/100);
#   LatEndMin=(lat_end-LatEndDeg*100)/60
#   LonEndMin=(lon_end-LonEndDeg*100)/60
# 
#   lat_start2= LatStartDeg + LatStartMin
#   lon_start2 = LonStartDeg + LonStartMin
#   lat_end2 = LatEndDeg + LatEndMin
#   lon_end2 = LonEndDeg + LonEndMin
#   Data$lat_start = lat_start2
#   Data$lon_start = lon_start2
#   Data$lat_end = lat_end2
#   Data$lon_end = lon_end2
# 
#   # quadrant 3
#   Data[Data$SHOOTING_QUADRANT == 3, "lat_start"] <- -1 * Data[Data$SHOOTING_QUADRANT == 3, "lat_start"]
#   Data[Data$HAULING_QUADRANT == 3, "lat_end"] <- -1 * Data[Data$HAULING_QUADRANT == 3, "lat_end"]
# 
#   # quadrant 5
#   Data[Data$SHOOTING_QUADRANT == 5, "lat_start"] <- -1 * Data[Data$SHOOTING_QUADRANT == 5, "lat_start"]
#   Data[Data$HAULING_QUADRANT == 5, "lat_end"] <- -1 * Data[Data$HAULING_QUADRANT == 5, "lat_end"]
# 
#   Data[Data$SHOOTING_QUADRANT == 5, "lon_start"] <- -1 * Data[Data$SHOOTING_QUADRANT == 5, "lon_start"]
#   Data[Data$HAULING_QUADRANT == 5, "lon_end"] <- -1 * Data[Data$HAULING_QUADRANT == 5, "lon_end"]
# 
#   # quadrant 7
#   Data[Data$SHOOTING_QUADRANT == 7, "lon_start"] <- -1 * Data[Data$SHOOTING_QUADRANT == 7, "lon_start"]
#   Data[Data$HAULING_QUADRANT == 7, "lon_end"] <- -1 * Data[Data$HAULING_QUADRANT == 7, "lon_end"]
# [MODIFICATO]:
  Data$lat_start = NA
  Data$lon_start = NA
  Data$lat_end = NA
  Data$lon_end = NA

  non_na_idx <- !is.na(Data$SHOOTING_LATITUDE) & !is.na(Data$SHOOTING_LONGITUDE) &
                !is.na(Data$HAULING_LATITUDE) & !is.na(Data$HAULING_LONGITUDE)

  if (any(non_na_idx)) {
    lat_start=Data$SHOOTING_LATITUDE[non_na_idx]
    lon_start= Data$SHOOTING_LONGITUDE[non_na_idx]
    lat_end=Data$HAULING_LATITUDE[non_na_idx]
    lon_end= Data$HAULING_LONGITUDE[non_na_idx]
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
    
    Data$lat_start[non_na_idx] = lat_start2
    Data$lon_start[non_na_idx] = lon_start2
    Data$lat_end[non_na_idx] = lat_end2
    Data$lon_end[non_na_idx] = lon_end2
  }

  # quadrant 3
  idx_3_shoot <- which(Data$SHOOTING_QUADRANT == 3 & !is.na(Data$lat_start))
  idx_3_haul  <- which(Data$HAULING_QUADRANT == 3 & !is.na(Data$lat_end))
  if (length(idx_3_shoot) > 0) {
    Data$lat_start[idx_3_shoot] <- -1 * Data$lat_start[idx_3_shoot]
  }
  if (length(idx_3_haul) > 0) {
    Data$lat_end[idx_3_haul] <- -1 * Data$lat_end[idx_3_haul]
  }

  # quadrant 5
  idx_5_shoot <- which(Data$SHOOTING_QUADRANT == 5 & !is.na(Data$lat_start))
  idx_5_haul  <- which(Data$HAULING_QUADRANT == 5 & !is.na(Data$lat_end))
  if (length(idx_5_shoot) > 0) {
    Data$lat_start[idx_5_shoot] <- -1 * Data$lat_start[idx_5_shoot]
    Data$lon_start[idx_5_shoot] <- -1 * Data$lon_start[idx_5_shoot]
  }
  if (length(idx_5_haul) > 0) {
    Data$lat_end[idx_5_haul] <- -1 * Data$lat_end[idx_5_haul]
    Data$lon_end[idx_5_haul] <- -1 * Data$lon_end[idx_5_haul]
  }

  # quadrant 7
  idx_7_shoot <- which(Data$SHOOTING_QUADRANT == 7 & !is.na(Data$lon_start))
  idx_7_haul  <- which(Data$HAULING_QUADRANT == 7 & !is.na(Data$lon_end))
  if (length(idx_7_shoot) > 0) {
    Data$lon_start[idx_7_shoot] <- -1 * Data$lon_start[idx_7_shoot]
  }
  if (length(idx_7_haul) > 0) {
    Data$lon_end[idx_7_haul] <- -1 * Data$lon_end[idx_7_haul]
  }

  return(Data)
}