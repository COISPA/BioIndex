#' Merge Ta-TB and TA-TC tables
#'
#' @param ta MEDITS or MEDITS-like TA table
#' @param tb MEDITS or MEDITS-like TB table
#' @param tc MEDITS or MEDITS-like TC table
#' @param species species rubin code (MEDITS format, e.g. "MERLMER")
#' @param country country code as reported in MEDITS format. "all" code to perform the analysis on all the countries of the same GSA
#' @param strata data frame of the stratification scheme adopet by the MEDITS survey
#' @param wd working directory
#' @param save boolean. If TRUE the plot is saved in the user defined working directory (wd)
#' @param verbose boolean. If TRUE a message is printed
#' @return A list of two data frames is returned. The first element contains the TA-TB merged tables, while the second element contains the TA-TC merged tables
#' @examples
#' m <- merge_TATBTC(ta=TA, tb=TB, tc=TC, species="MERLMER", country="all", wd=tempdir(), verbose=TRUE)
#' mTATB <- m[[1]]
#' mTATC <- m[[2]]
#'
#' @export merge_TATBTC
#' @importFrom hms hms
#' @importFrom utils write.table

merge_TATBTC <- function(ta, tb, tc, species, country="all", strata=BioIndex::strata_scheme, wd=NA, save=TRUE, verbose=TRUE) {

  if (FALSE) {
    # library(BioIndex)
    verbose=TRUE
    wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\R_BioIndex_3.3_(in update)"
    ta <- read.table(paste(wd,"input/TA/TA GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
    tb <- read.table(paste(wd,"input/TB/TB GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
    tc <- read.table(paste(wd,"input/TC/TC GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
    country="all"

    # load("D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\- Pacchetto R BioIndex -\\BioIndex\\data\\TA_cols.rda")
    # load("D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\- Pacchetto R BioIndex -\\BioIndex\\data\\TB_cols.rda")
    # load("D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\- Pacchetto R BioIndex -\\BioIndex\\data\\TC_cols.rda")
    species <- "ARISFOL"




    ta <- read.table(file.path(wd,"input","medits_ta_1.csv"),sep=";",header=TRUE)
    tb <- read.table(file.path(wd,"input","medits_tb_1.csv"),sep=";",header=TRUE)
    tc <- read.table(file.path(wd,"input","medits_tc_1.csv"),sep=";",header=TRUE)


    m <- merge_TATBTC(ta, tb, tc, species="PAPELON", country="all", wd=wd, verbose=TRUE)
  }


  if (is.na(wd) & save) {
    save =FALSE
    if (verbose){
      message("Missing working directory. Results are not saved in the local folder.")
    }
  }

  strata_scheme <- strata

  check_date_haul <- check_dictionary <-
    check_hauls_TBTA <- check_numeric_range <- convert_coordinates <- hms <-
    write.table <- NULL

  TA_cols <- BioIndex::TA_cols
  TB_cols <- BioIndex::TB_cols
  TC_cols <- BioIndex::TC_cols



  species[2] <- substr(species,5,7)
  species[1] <- substr(species[1],1,4)
  sspp <- paste(species[1],species[2], sep="")
  GENERE <- species[1]
  SPECIE <- species[2]

  # if (verbose){
  #     message("- metaDB preparation -\n")
  # }
  TA <- ta
  TB <- tb
  TC <- tc


  #########################
  ###   filtro COUNTRY  ###
  #########################
  TA_country <- sort(unique(as.character(TA$COUNTRY)))
  TB_country <- sort(unique(as.character(TB$COUNTRY)))
  TC_country <- sort(unique(as.character(TC$COUNTRY)))
  l_country <- length(TA_country)

  for (i in 1:l_country){
    if (!(TA_country[i] %in% TB_country)){warning(paste("The country ",TA_country[i]," is not present in the TB file.", sep=""))}
    if (!(TA_country[i] %in% TC_country)){warning(paste("The country ",TA_country[i]," is not present in the TC file.", sep=""))}
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


  TA <- TA[TA$COUNTRY %in% country_analysis , ] #  & TA$VALIDITY=="V"
  TB <- TB[TB$COUNTRY %in% country_analysis , ]
  TC <- TC[TC$COUNTRY %in% country_analysis , ]


  id_TA <- data.frame(id = paste(TA$AREA,TA$COUNTRY,TA$YEAR,"_", TA$VESSEL, TA$MONTH, TA$DAY,"_", TA$HAUL_NUMBER, sep = ""))
  id_TB <- data.frame(id = paste(TB$AREA,TB$COUNTRY,TB$YEAR,"_", TB$VESSEL, TB$MONTH, TB$DAY,"_", TB$HAUL_NUMBER, sep = ""))
  id_TC <- data.frame(id = paste(TC$AREA,TC$COUNTRY,TC$YEAR,"_", TC$VESSEL, TC$MONTH, TC$DAY,"_", TC$HAUL_NUMBER, sep = ""))



  id_invalid <- id_TA$id[which(TA$VALIDITY == "I")]

  ##----------------
  ## RoMEBS checks
  ##----------------
  suffix <- paste(as.character(Sys.Date()),format(Sys.time(), "_time_h%Hm%Ms%OS0"),sep="")
  yyy <- sort(unique(TA$YEAR))
  yy <- 28
  for (yy in 1:length(yyy)) {
    # print(yy)
    TA_year <- TA[TA$YEAR==yyy[yy], ]
    TB_year <- TB[TB$YEAR==yyy[yy], ]
    TC_year <- TC[TC$YEAR==yyy[yy], ]

    ## TA checks
    Field = "SHOOTING_TIME"
    Values = seq(0,2400,1)
    if (!(check_dictionary(ResultData = TA_year, Field, Values, year=yyy[yy], wd=paste(wd, "/output/", sep=""), suffix=suffix))) {
      stop("SHOOTING_TIME value out of the allowed range. Please check logfile in ~/output/Logliles")
    }
    Field = "HAULING_TIME"
    Values = seq(0,2400,1)
    if (!(check_dictionary(ResultData = TA_year, Field, Values, year=yyy[yy], wd=paste(wd, "/output/", sep=""), suffix=suffix))) {
      stop("HAULING_TIME value out of the allowed range. Please check logfile in ~/output/Logliles")
    }
    Field = "WING_OPENING"
    Values = c(30, 50:250)
    if (!(check_dictionary(ResultData = TA_year, Field, Values, year=yyy[yy], wd=paste(wd, "/output/", sep=""), suffix=suffix) ) ) {
      stop("WING_OPENING value out of the allowed range. Please check logfile in ~/output/Logliles")
    }
    Field = "SHOOTING_LATITUDE"
    Values = c(3020,4730)
    if (!(check_numeric_range(Data = TA_year, Field, Values, year=yyy[yy], wd=paste(wd, "/output/", sep=""), suffix=suffix))) {
      stop("SHOOTING_LATITUDE in TA out of the allowed range. Please check logfile in ~/output/Logliles")
    }
    Field = "HAULING_LATITUDE"
    Values = c(3020,4730)
    if (!(check_numeric_range(Data = TA_year, Field, Values, year=yyy[yy], wd=paste(wd, "/output/", sep=""), suffix=suffix))) {
      stop("HAULING_LATITUDE in TA out of the allowed range. Please check logfile in ~/output/Logliles")
    }
    Field = "SHOOTING_LONGITUDE"
    Values = c(0,4200.00)
    if (!(check_numeric_range(Data = TA_year, Field, Values, year=yyy[yy], wd=paste(wd, "/output/", sep=""), suffix=suffix))) {
      stop("SHOOTING_LONGITUDE in TA out of the allowed range. Please check logfile in ~/output/Logliles")
    }
    Field = "HAULING_LONGITUDE"
    Values = c(0,4200.00)
    if (!(check_numeric_range(Data = TA_year, Field, Values, year=yyy[yy], wd=paste(wd, "/output/", sep=""), suffix=suffix))) {
      stop("HAULING_LONGITUDE in TA out of the allowed range. Please check logfile in ~/output/Logliles")
    }
    ## RoMEBS check of Date consistency between TA and TB
    if ( !(check_date_haul(TA_year,TB_year, year=yyy[yy],wd=paste(wd, "/output/", sep=""), suffix=suffix)) ) {
      stop("Date reported in TB not consistent with the date reported in TA. Please check logfile in ~/output/Logliles")
    }

    if ( !(check_hauls_TBTA(TA_year,TB_year, year=yyy[yy],wd=paste(wd, "/output/", sep=""), suffix=suffix)) ) {
      stop("Haul in TB not reported in TA. Please check logfile in ~/output/Logliles")
    }
    ## RoMEBS check of Date consistency between TA and TC
    if ( !(check_date_haul(TA_year,TC_year, year=yyy[yy],wd=paste(wd, "/output/", sep=""), suffix=suffix)) ) {
      stop("Date reported in TC not consistent with the date reported in TA. Please check logfile in ~/output/Logliles")
    }
  }


  colnames(TA)[which(colnames(TA) == "AREA")] <- "GSA"
  colnames(TB)[which(colnames(TB) == "AREA")] <- "GSA"
  colnames(TC)[which(colnames(TC) == "AREA")] <- "GSA"
  TA_merge <- cbind(id_TA,TA)
  TB_merge <- cbind(id_TB,TB)
  TB_merge$GENUS <- as.character(TB_merge$GENUS)
  TB_merge$SPECIES <- as.character(TB_merge$SPECIES)
  TC_merge <- cbind(id_TC,TC)
  TC_merge$GENUS <- as.character(TC_merge$GENUS)
  TC_merge$SPECIES <- as.character(TC_merge$SPECIES)

  #remove invalid hauls
  TA_merge <- TA_merge[!TA_merge$id %in% id_invalid, ]
  TB_merge <- TB_merge[!TB_merge$id %in% id_invalid, ]
  TC_merge <- TC_merge[!TC_merge$id %in% id_invalid, ]

  TA_merge <- TA_merge[,which(colnames(TA_merge) %in% TA_cols)]
  TB_merge <- TB_merge[TB_merge$GENUS == species[1] & TB_merge$SPECIES == species[2], which(colnames(TB_merge) %in% TB_cols)]
  TC_merge <- TC_merge[TC_merge$GENUS == species[1] & TC_merge$SPECIES == species[2], which(colnames(TC_merge) %in% TC_cols)]

  ########################
  ##    MERGE TA- TB    ##
  ########################
  if(verbose){
    message("- Merging TA-TB files")
  }
  merge_TATB <- merge(TA_merge, TB_merge, by.x = "id", by.y = "id", all.x = T)
  merge_TATB$MEDITS_CODE <- as.character(paste(merge_TATB$GENUS, merge_TATB$SPECIES))

  l_TATB <- length (merge_TATB[,1])
  i=1
  for (i in 1:l_TATB){
    if (merge_TATB[i, "MEDITS_CODE"] == "NA NA"){
      merge_TATB[i, "MEDITS_CODE"              ] <- "NA"
      merge_TATB[i, "GENUS"                    ] <- -1
      merge_TATB[i, "SPECIES"                  ] <- -1
      merge_TATB[i, "TOTAL_WEIGHT_IN_THE_HAUL" ] <- 0
      merge_TATB[i, "TOTAL_NUMBER_IN_THE_HAUL" ] <- 0
      merge_TATB[i, "NB_OF_FEMALES"            ] <- 0
      merge_TATB[i, "NB_OF_MALES"              ] <- 0
      merge_TATB[i, "NB_OF_UNDETERMINED"       ] <- 0
    }
  }

  Data <- merge_TATB
  coord <- convert_coordinates(merge_TATB)

  merge_TATB$MEAN_LATITUDE <- (merge_TATB$SHOOTING_LATITUDE+merge_TATB$HAULING_LATITUDE)/2
  merge_TATB$MEAN_LATITUDE_DEC <- (coord$lat_start+coord$lat_end)/2
  merge_TATB$MEAN_LONGITUDE <- (merge_TATB$SHOOTING_LONGITUDE+merge_TATB$HAULING_LONGITUDE)/2
  merge_TATB$MEAN_LONGITUDE_DEC <- (coord$lon_start+coord$lon_end)/2
  merge_TATB$MEAN_DEPTH <- (merge_TATB$SHOOTING_DEPTH+merge_TATB$HAULING_DEPTH)/2
  merge_TATB$SWEPT_AREA <- merge_TATB$DISTANCE * merge_TATB$WING_OPENING/10000000
  i=1
  hour_shooting <- 0
  min_shooting <- 0
  hour_hauling <- 0
  min_hauling <- 0

  strata_scheme <- strata_scheme[strata_scheme$GSA == unique(merge_TATB$GSA) & strata_scheme$COUNTRY %in% as.character(unique(merge_TATB$COUNTRY)), ]

i=1
  for (i in 1:l_TATB){
j=1
    for (j in 1:length(strata_scheme$CODE)){

      if(strata_scheme[j,3]==1){
        if  (floor(merge_TATB$MEAN_DEPTH[i]) >= strata_scheme[j,4] & floor(merge_TATB$MEAN_DEPTH[i]) <= strata_scheme[j,5]) {merge_TATB$STRATUM_CODE[i] <- strata_scheme[j,3]}
      } else {
        if  (floor(merge_TATB$MEAN_DEPTH[i]) > strata_scheme[j,4] & floor(merge_TATB$MEAN_DEPTH[i]) <= strata_scheme[j,5]) {merge_TATB$STRATUM_CODE[i] <- strata_scheme[j,3]}
      }

    }

    if (nchar(merge_TATB$SHOOTING_TIME[i])==4) {hour_shooting[i] <- substr(merge_TATB$SHOOTING_TIME[i],1,2); min_shooting[i] <- substr(merge_TATB$SHOOTING_TIME[i],3,4)} else {
      if (nchar(merge_TATB$SHOOTING_TIME[i])==3) {hour_shooting[i] <- substr(merge_TATB$SHOOTING_TIME[i],1,1); min_shooting[i] <- substr(merge_TATB$SHOOTING_TIME[i],2,3)}}

    if (nchar(merge_TATB$HAULING_TIME[i])==4) {hour_hauling[i] <- substr(merge_TATB$HAULING_TIME[i],1,2); min_hauling[i] <- substr(merge_TATB$HAULING_TIME[i],3,4)} else {
      if (nchar(merge_TATB$HAULING_TIME[i])==3) {hour_hauling[i] <- substr(merge_TATB$HAULING_TIME[i],1,1); min_hauling[i] <- substr(merge_TATB$HAULING_TIME[i],2,3)}}
  }

  hms_shooting <- hms( rep (0, length(hour_shooting)), as.numeric(min_shooting), as.numeric(hour_shooting))
  hms_hauling  <- hms( rep (0, length(hour_shooting)), as.numeric(min_hauling), as.numeric(hour_hauling))
  duration <-  as.numeric(hms_hauling - hms_shooting)/3600

  merge_TATB$SHOOTING_TIME <- hms_shooting
  merge_TATB$HAULING_TIME  <- hms_hauling

  merge_TATB$N_h <- merge_TATB$TOTAL_NUMBER_IN_THE_HAUL/duration
  merge_TATB$N_km2 <- merge_TATB$TOTAL_NUMBER_IN_THE_HAUL/merge_TATB$SWEPT_AREA
  merge_TATB$kg_h <- merge_TATB$TOTAL_WEIGHT_IN_THE_HAUL/duration/1000
  merge_TATB$kg_km2 <- merge_TATB$TOTAL_WEIGHT_IN_THE_HAUL/merge_TATB$SWEPT_AREA/1000

  merge_TATB$MEDITS_CODE <- as.character(merge_TATB$MEDITS_CODE)

  if (save){
    write.table(merge_TATB, paste(wd, "/output/mergeTATB_",species[1], species[2],".csv", sep=""), sep=";", row.names=F)
  }

  if(verbose){
  message("TA-TB files correctly merged")
  message(paste("Merge TA-TB files saved in the following folder: '",wd, "/output/mergeTATB_",species[1], species[2],".csv'\n", sep=""))
}
  ########################
  ##    MERGE TA- TC    ##
  ########################
  if(verbose){
      message("- Merging TA-TC files")
  }

  merge_TATC <- merge(TA_merge, TC_merge, by.x = "id", by.y = "id", all.x = T, all.y = T)
  merge_TATC$MEDITS_CODE <- paste(merge_TATC$GENUS, merge_TATC$SPECIES)

  merge_TATC$LENGTH_CLASSES_CODE <- as.character(merge_TATC$LENGTH_CLASSES_CODE)
  merge_TATC$SEX <- as.character(merge_TATC$SEX)
  merge_TATC$MATURITY <- as.character(merge_TATC$MATURITY)
  merge_TATC$MATSUB <- toupper(as.character(merge_TATC$MATSUB))
  l_TATC <- length (merge_TATC[,1])
  i=2
  for (i in 1:l_TATC){
    if (merge_TATC[i, "MEDITS_CODE"] == "NA NA"){
      merge_TATC[i, "MEDITS_CODE"              ] <- -1
      merge_TATC[i, "GENUS"                    ] <- -1
      merge_TATC[i, "SPECIES"                  ] <- -1
      merge_TATC[i, "SEX"                      ] <- -1
      merge_TATC[i, "LENGTH_CLASSES_CODE"      ] <- -1
      merge_TATC[i, "WEIGHT_OF_THE_FRACTION"   ] <-  0
      merge_TATC[i, "WEIGHT_OF_THE_SAMPLE_MEASURED"   ] <- 0
      merge_TATC[i, "MATURITY"                 ] <- -1
      merge_TATC[i, "MATSUB"                   ] <- -1
      merge_TATC[i, "LENGTH_CLASS"             ] <- -1
      merge_TATC[i, "NUMBER_OF_INDIVIDUALS_IN_THE_LENGTH_CLASS_AND_MATURITY_STAGE"] <- 0
    }

    merge_TATC[i,]

  }

  coord <- convert_coordinates(merge_TATC)

  merge_TATC$MEAN_LATITUDE <- (merge_TATC$SHOOTING_LATITUDE+merge_TATC$HAULING_LATITUDE)/2
  merge_TATC$MEAN_LATITUDE_DEC <- (coord$lat_start+coord$lat_end)/2
  merge_TATC$MEAN_LONGITUDE <- (merge_TATC$SHOOTING_LONGITUDE+merge_TATC$HAULING_LONGITUDE)/2
  merge_TATC$MEAN_LONGITUDE_DEC <- (coord$lon_start+coord$lon_end)/2
  merge_TATC$MEAN_DEPTH <- (merge_TATC$SHOOTING_DEPTH+merge_TATC$HAULING_DEPTH)/2
  merge_TATC$SWEPT_AREA <- merge_TATC$DISTANCE * merge_TATC$WING_OPENING/10000000

strata_scheme <- strata_scheme[strata_scheme$GSA == unique(merge_TATB$GSA) & strata_scheme$COUNTRY %in% as.character(unique(merge_TATB$COUNTRY)), ]


  i=1
  hour_shooting <- 0
  min_shooting <- 0
  hour_hauling <- 0
  min_hauling <- 0

  for (i in 1:l_TATC){

    for (j in 1:length(strata_scheme$CODE)){
      if  (merge_TATC$MEAN_DEPTH[i] > strata_scheme[j,4] & merge_TATC$MEAN_DEPTH[i] <= strata_scheme[j,5]) {merge_TATC$STRATUM_CODE[i] <- strata_scheme[j,3]}
    }

    if (nchar(merge_TATC$SHOOTING_TIME[i])==4) {hour_shooting[i] <- substr(merge_TATC$SHOOTING_TIME[i],1,2); min_shooting[i] <- substr(merge_TATC$SHOOTING_TIME[i],3,4)} else {
      if (nchar(merge_TATC$SHOOTING_TIME[i])==3) {hour_shooting[i] <- substr(merge_TATC$SHOOTING_TIME[i],1,1); min_shooting[i] <- substr(merge_TATC$SHOOTING_TIME[i],2,3)}}

    if (nchar(merge_TATC$HAULING_TIME[i])==4) {hour_hauling[i] <- substr(merge_TATC$HAULING_TIME[i],1,2); min_hauling[i] <- substr(merge_TATC$HAULING_TIME[i],3,4)} else {
      if (nchar(merge_TATC$HAULING_TIME[i])==3) {hour_hauling[i] <- substr(merge_TATC$HAULING_TIME[i],1,1); min_hauling[i] <- substr(merge_TATC$HAULING_TIME[i],2,3)}}
  }



  hms_shooting <- hms( rep (0, length(hour_shooting)), as.numeric(min_shooting), as.numeric(hour_shooting))
  hms_hauling  <- hms( rep (0, length(hour_shooting)), as.numeric(min_hauling), as.numeric(hour_hauling))
  duration <-  as.numeric(hms_hauling - hms_shooting)/3600

  merge_TATC$SHOOTING_TIME <- hms_shooting
  merge_TATC$HAULING_TIME  <- hms_hauling

  i=2
  for (i in 1:l_TATC){
    if(merge_TATC[i, "MATURITY"] == -1) {merge_TATC[i, "MATURITY_STAGE"] <- -1 }
    if(merge_TATC[i, "MATURITY"] == 0) {merge_TATC[i, "MATURITY_STAGE"] <- 0 }
    if(merge_TATC[i, "MATURITY"] == 1) {merge_TATC[i, "MATURITY_STAGE"] <- 1 }
    if(merge_TATC[i, "MATURITY"] == 2 & !(merge_TATC[i, "MATSUB"] %in% c("A","B","C","D","E","F"))) {merge_TATC[i, "MATURITY_STAGE"] <- 2 } else {
      if (merge_TATC[i, "MATURITY"] == 2 & merge_TATC[i, "MATSUB"] %in% c("A","B","C","D","E","F")) {
        merge_TATC[i, "MATURITY_STAGE"] <- paste(merge_TATC[i, "MATURITY"],merge_TATC[i, "MATSUB"],sep="")}
    }
    if(merge_TATC[i, "MATURITY"] == 3 & !(merge_TATC[i, "MATSUB"] %in% c("A","B","C","D","E","F"))) {merge_TATC[i, "MATURITY_STAGE"] <- 3 } else {
      if (merge_TATC[i, "MATURITY"] == 3 & merge_TATC[i, "MATSUB"] %in% c("A","B","C","D","E","F")) {
        merge_TATC[i, "MATURITY_STAGE"] <- paste(merge_TATC[i, "MATURITY"],merge_TATC[i, "MATSUB"],sep="")}
    }
    if(merge_TATC[i, "MATURITY"] == 4 & !(merge_TATC[i, "MATSUB"] %in% c("A","B","C","D","E","F"))) {merge_TATC[i, "MATURITY_STAGE"] <- 4 } else {
      if (merge_TATC[i, "MATURITY"] == 4 & merge_TATC[i, "MATSUB"] %in% c("A","B","C","D","E","F")) {
        merge_TATC[i, "MATURITY_STAGE"] <- paste(merge_TATC[i, "MATURITY"],merge_TATC[i, "MATSUB"],sep="")}
    }

    if (merge_TATC[i, "WEIGHT_OF_THE_FRACTION"] > 0){merge_TATC[i, "RAISING_FACTOR"] <- merge_TATC[i, "WEIGHT_OF_THE_FRACTION"] / merge_TATC[i, "WEIGHT_OF_THE_SAMPLE_MEASURED"]
    } else {merge_TATC[i, "RAISING_FACTOR"] <- 0}
    if (merge_TATC[i, "RAISING_FACTOR"] >0 & merge_TATC[i, "RAISING_FACTOR"] < 1) {merge_TATC[i, "RAISING_FACTOR"] <- 1}
  }

  merge_TATC$N_h <- merge_TATC$NUMBER_OF_INDIVIDUALS_IN_THE_LENGTH_CLASS_AND_MATURITY_STAGE*merge_TATC$RAISING_FACTOR/duration
  merge_TATC$N_km2 <- merge_TATC$NUMBER_OF_INDIVIDUALS_IN_THE_LENGTH_CLASS_AND_MATURITY_STAGE*merge_TATC$RAISING_FACTOR/merge_TATC$SWEPT_AREA
  merge_TATC$kg_h <- merge_TATC$WEIGHT_OF_THE_SAMPLE_MEASURED/duration/1000
  merge_TATC$kg_km2 <- merge_TATC$WEIGHT_OF_THE_SAMPLE_MEASURED/merge_TATC$SWEPT_AREA/1000

  if (save){
    write.table(merge_TATC, paste(wd, "/output/mergeTATC_",species[1], species[2],".csv", sep=""), sep=";", row.names=F)
  }

  if(verbose){
    message("TA-TC files correctly merged")
    message(paste("Merge TA-TC files saved in the following folder: '",wd, "/output/mergeTATC_",species[1], species[2],".csv'\n", sep=""))
  }

return(list(merge_TA_TB = merge_TATB, merge_TA_TC = merge_TATC))
}
