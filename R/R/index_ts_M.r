#' Estimation of abundance indices for males
#'
#' @param mTATB data frame
#' @param GSA reference GSA for the analysis
#' @param country_analysis vector of reference countries for the analysis
#' @param depth_range range of depth strata to perform the analysis (min, max)
#' @param strata_scheme data frame of the stratification scheme
#' @param stratification data frame of strata surface area
#' @param wd working directory
#' @param save boolean. If TRUE the plot is saved in the user defined working directory (wd)
#' @importFrom gridExtra grid.arrange
#' @importFrom reshape2 melt
#' @importFrom stringr str_split
#' @importFrom grDevices dev.off tiff
#' @importFrom graphics legend lines par
#' @export
index_ts_M <- function(mTATB, GSA, country_analysis, depth_range, strata_scheme, stratification, wd=NA, save=TRUE) {

    MEAN_DEPTH <- strata <- value <- year <- NULL
    merge_TATB <- mTATB

    merge_TATB$N_km2 <- merge_TATB$NB_OF_MALES / merge_TATB$SWEPT_AREA

    #-----------------------
    # Abundance INDICES
    #-----------------------

    #------------------------
    # Variable definitions #
    dependent <- "abundance"
    dep_text <-expression(paste("Abundance ", (n/km^2), sep=" "))
    varcol <- which(colnames(merge_TATB)=="N_km2")
    col_response <- which(colnames(merge_TATB)=="N_km2")
    positive_hauls<- 0
    total_hauls<- 0

    GENERE <- as.character(unique(merge_TATB$GENUS)[unique(merge_TATB$GENUS) != -1])
    SPECIE <- as.character(unique(merge_TATB$SPECIES)[unique(merge_TATB$SPECIES) != -1])
    sspp <- paste(GENERE,SPECIE, sep="")
    GSA <- unique(merge_TATB$GSA)
    species <- sspp
    #------------------------


    #------------------------------------------
    # selection of depth range
    data <- merge_TATB
    depth <- data.frame(matrix(NA,ncol=4, nrow=length(strata_scheme$CODE)))
    colnames(depth) <- c("strata", "min", "max", "bul")
    depth$strata <- strata_scheme$CODE # c(1,2,3,4,5)
    depth$min <-strata_scheme$MIN_DEPTH  # c(10,50,100,200,500)
    depth$max <-strata_scheme$MAX_DEPTH  # c(50,100,200,500,800)
    # print("Select the depth range for the analysis")
    # depth_range <- dlgInput("Depth range for the analysis: ", default="10,100",Sys.info()[""])$res
    # depth_range <- data.frame(range = strsplit(depth_range, ",")); depth_range <- as.numeric(as.character(depth_range[,1]))
    if (depth_range[2] != 800) {depth_range[2] <- depth_range[2]}
    data <-  data[data$MEAN_DEPTH>depth_range[1] & data$MEAN_DEPTH<=depth_range[2],]
    strata_range <- seq(which(depth[,"min"]==depth_range[1]),which(depth[,"max"]==depth_range[2]),1)
    depth[depth$strata %in% strata_range, "bul"] <- T
    depth[!(depth$strata %in% strata_range), "bul"] <- F

    analysis_stratum1 <- NA
    analysis_stratum2 <- NA
    analysis_stratum3 <- NA
    analysis_stratum4 <- NA
    analysis_stratum5 <- NA
    analysis_stratum6 <- NA

    analysis_stratum1 <- depth[1,4]
    analysis_stratum2 <- depth[2,4]
    analysis_stratum3 <- depth[3,4]
    analysis_stratum4 <- depth[4,4]
    analysis_stratum5 <- depth[5,4]
    analysis_stratum6 <- depth[6,4]

    if (is.na(analysis_stratum1)){depth[1,4]<- F; analysis_stratum1 <- F}
    if (is.na(analysis_stratum2)){depth[2,4]<- F; analysis_stratum2 <- F}
    if (is.na(analysis_stratum3)){depth[3,4]<- F; analysis_stratum3 <- F}
    if (is.na(analysis_stratum4)){depth[4,4]<- F; analysis_stratum4 <- F}
    if (is.na(analysis_stratum5)){depth[5,4]<- F; analysis_stratum5 <- F}
    if (is.na(analysis_stratum6)){depth[6,4]<- F; analysis_stratum6 <- F}

    #------------------------------------------


    ddd <- merge_TATB
    year_range <- data.frame(year = sort(unique(merge_TATB$YEAR)))

    res_table_cala2 <- data.frame(year=sort(unique(merge_TATB$YEAR)))
    res_table_cala2$index_1 <- NA
    res_table_cala2$index_2 <- NA
    res_table_cala2$index_3 <- NA
    res_table_cala2$index_4 <- NA
    res_table_cala2$index_5 <- NA
    res_table_cala2$index_6 <- NA

    se_table2  <-  data.frame(year=sort(unique(merge_TATB$YEAR)))
    se_table2$sum_1 <- NA
    se_table2$sum_2 <- NA
    se_table2$sum_3 <- NA
    se_table2$sum_4 <- NA
    se_table2$sum_5 <- NA
    se_table2$sum_6 <- NA
    se_table2$sd <- NA

    # stratification <- read.table(paste(wd,"/scripts/utilities/stratification scheme.csv", sep=""), sep=";", header=T)
    stratification$SURF <- as.numeric(stratification$SURF)

    if (analysis_stratum1 == T){area_s1 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==1 & stratification$COUNTRY %in% country_analysis,5])} else {area_s1 <- 0}
    if (analysis_stratum2 == T){area_s2 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==2 & stratification$COUNTRY %in% country_analysis,5])} else {area_s2 <- 0}
    if (analysis_stratum3 == T){area_s3 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==3 & stratification$COUNTRY %in% country_analysis,5])} else {area_s3 <- 0}
    if (analysis_stratum4 == T){area_s4 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==4 & stratification$COUNTRY %in% country_analysis,5])} else {area_s4 <- 0}
    if (analysis_stratum5 == T){area_s5 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==5 & stratification$COUNTRY %in% country_analysis,5])} else {area_s5 <- 0}
    if (analysis_stratum6 == T){area_s6 <- sum(stratification[stratification$GSA == GSA & stratification$CODE==6 & stratification$COUNTRY %in% country_analysis,5])} else {area_s6 <- 0}


    peso_s1 <- area_s1/sum(area_s1,area_s2,area_s3,area_s4,area_s5,area_s6)
    peso_s2 <- area_s2/sum(area_s1,area_s2,area_s3,area_s4,area_s5,area_s6)
    peso_s3 <- area_s3/sum(area_s1,area_s2,area_s3,area_s4,area_s5,area_s6)
    peso_s4 <- area_s4/sum(area_s1,area_s2,area_s3,area_s4,area_s5,area_s6)
    peso_s5 <- area_s5/sum(area_s1,area_s2,area_s3,area_s4,area_s5,area_s6)
    peso_s6 <- area_s6/sum(area_s1,area_s2,area_s3,area_s4,area_s5,area_s6)

    i=1
    for(i in 1:length(res_table_cala2[,1])) {

        data2 <- ddd[ddd$YEAR == res_table_cala2[i,1] , ]

        if (analysis_stratum1 == T) {s1 <- data2[floor(data2$MEAN_DEPTH) >= depth[1,2] & floor(data2$MEAN_DEPTH) <= depth[1,3], ]}
        if (analysis_stratum2 == T) {s2 <- data2[floor(data2$MEAN_DEPTH) > depth[2,2] & floor(data2$MEAN_DEPTH) <= depth[2,3], ]}
        if (analysis_stratum3 == T) {s3 <- data2[floor(data2$MEAN_DEPTH) > depth[3,2] & floor(data2$MEAN_DEPTH) <= depth[3,3], ]}
        if (analysis_stratum4 == T) {s4 <- data2[floor(data2$MEAN_DEPTH) > depth[4,2] & floor(data2$MEAN_DEPTH) <= depth[4,3], ]}
        if (analysis_stratum5 == T) {s5 <- data2[floor(data2$MEAN_DEPTH) > depth[5,2] & floor(data2$MEAN_DEPTH) <= depth[5,3], ]}
        if (analysis_stratum6 == T) {s6 <- data2[floor(data2$MEAN_DEPTH) > depth[6,2] & floor(data2$MEAN_DEPTH) <= depth[6,3], ]}





        if (analysis_stratum1 == T){
            s_n <- sum(s1[!is.na(s1$NB_OF_MALES),"NB_OF_MALES"])
            s_a <- sum(s1$SWEPT_AREA)
            res_table_cala2[i,2] <- s_n/s_a

            mean_ind <- mean(s1[!is.na(s1$N_km2), "N_km2" ])
            sq <- (s1[!is.na(s1$N_km2),"N_km2"] - mean_ind)^2
            sqA <- sq*s1$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s1[,1])-1)
            f1 <- s_a/area_s1
            se_table2[i,2] <- (((peso_s1)^2 * sum_sqA)/ s_a)* (1-f1)

        } else {area_s1 = 0}

        if (analysis_stratum2 == T){
            s_n <- sum(s2[!is.na(s2$NB_OF_MALES),"NB_OF_MALES"])
            s_a <- sum(s2$SWEPT_AREA)
            res_table_cala2[i,3] <- s_n/s_a

            mean_ind <- mean(s2[!is.na(s2$N_km2), "N_km2" ])
            sq <- (s2[!is.na(s2$N_km2),"N_km2"] - mean_ind)^2
            sqA <- sq*s2$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s2[,1])-1)
            f2 <- s_a/area_s2
            se_table2[i,3] <- (((peso_s2)^2 * sum_sqA)/ s_a)* (1-f2)
        } else {area_s2 = 0}

        if (analysis_stratum3 == T){
            s_n <- sum(s3[!is.na(s3$NB_OF_MALES),"NB_OF_MALES"])
            s_a <- sum(s3$SWEPT_AREA)
            res_table_cala2[i,4] <- s_n/s_a

            mean_ind <- mean(s3[!is.na(s3$N_km2), "N_km2" ])
            sq <- (s3[!is.na(s3$N_km2),"N_km2"] - mean_ind)^2
            sqA <- sq*s3$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s3[,1])-1)
            f3 <- s_a/area_s3
            se_table2[i,4] <- (((peso_s3)^2 * sum_sqA)/ s_a)* (1-f3)
        } else {area_s3 = 0}

        if (analysis_stratum4 == T){
            s_n <- sum(s4[!is.na(s4$NB_OF_MALES),"NB_OF_MALES"])
            s_a <- sum(s4$SWEPT_AREA)
            res_table_cala2[i,5] <- s_n/s_a

            mean_ind <- mean(s4[!is.na(s4$N_km2), "N_km2" ])
            sq <- (s4[!is.na(s4$N_km2),"N_km2"] - mean_ind)^2
            sqA <- sq*s4$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s4[,1])-1)
            f4 <- s_a/area_s4
            se_table2[i,5] <- (((peso_s4)^2 * sum_sqA)/ s_a)* (1-f4)
        } else {area_s4 = 0}

        if (analysis_stratum5 == T){
            #index computation
            s_n <- sum(s5[!is.na(s5$NB_OF_MALES),"NB_OF_MALES"])
            s_a <- sum(s5$SWEPT_AREA)
            res_table_cala2[i,6] <- s_n/s_a

            #se computation
            mean_ind <- mean(s5[!is.na(s5$N_km2), "N_km2" ])
            sq <- (s5[!is.na(s5$N_km2),"N_km2"] - mean_ind)^2
            sqA <- sq*s5$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s5[,1])-1)
            f5 <- s_a/area_s5
            se_table2[i,6] <- (((peso_s5)^2 * sum_sqA)/ s_a)* (1-f5)
        } else {area_s5 = 0}

        if (analysis_stratum6 == T){
            #index computation
            s_n <- sum(s6[!is.na(s6$NB_OF_MALES),"NB_OF_MALES"])
            s_a <- sum(s6$SWEPT_AREA)
            res_table_cala2[i,7] <- s_n/s_a

            #se computation
            mean_ind <- mean(s6[!is.na(s6$N_km2), "N_km2" ])
            sq <- (s6[!is.na(s6$N_km2),"N_km2"] - mean_ind)^2
            sqA <- sq*s6$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s6[,1])-1)
            f6 <- s_a/area_s6
            se_table2[i,7] <- (((peso_s6)^2 * sum_sqA)/ s_a)* (1-f6)
        } else {area_s6 = 0}

        # source(paste(wd, "/scripts/Index_estimation_cala.R", sep=""), encoding = 'UTF-8')
        positive_hauls[i] <- length(merge_TATB[merge_TATB$NB_OF_MALES>0 & merge_TATB$YEAR == res_table_cala2[i,1],"NB_OF_MALES"])
        total_hauls[i] <- length(merge_TATB[merge_TATB$YEAR == res_table_cala2[i,1],"HAUL_NUMBER"])

        sum_res_est <- c(res_table_cala2[i,2]*peso_s1,res_table_cala2[i,3]*peso_s2,res_table_cala2[i,4]*peso_s3,res_table_cala2[i,5]*peso_s4,res_table_cala2[i,6]*peso_s5,res_table_cala2[i,7]*peso_s6)
        res_table_cala2[i, 8]<- sum(sum_res_est[!is.na(sum_res_est)])# /sum(area_s1,area_s2,area_s3,area_s4,area_s5)
        colnames(res_table_cala2) <- c("year", "stratum 1","stratum 2", "stratum 3", "stratum 4", "stratum 5", "stratum 6", "Indices")
        se_table2[i, "sd"] <-sqrt(rowSums(se_table2[i, 2:7], na.rm = T))
    }

    timeseries <- data.frame(year = res_table_cala2[,1], abundance = res_table_cala2[,8], sd= se_table2[, "sd"])
    timeseries$CV <- timeseries[,3] / timeseries[,2]
    timeseries$invCV <- 1/timeseries$CV
    timeseries$positive_hauls_perc <- positive_hauls / total_hauls * 100
    colnames(res_table_cala2) <- c("year", "stratum 1","stratum 2", "stratum 3", "stratum 4", "stratum 5", "stratum 6", "Indices")

    # plot timeseries of mean indices
    main <- paste(sspp,"_GSA",GSA,"_(",dependent,")-MALES-Random_Stratified_Sampling_",depth_range[1],"-",depth_range[2], " m", sep="")
    main.lab <- paste(sspp," GSA",GSA," (",dependent,")-MALES-RSS ",depth_range[1],"-",depth_range[2], " m", sep="")
    max_index <- max(timeseries[,2]) + (max(timeseries[!is.na(timeseries$sd),"sd"])*1.2)
    if(save){
        write.table(timeseries, paste(wd,"/output/",main,"_Timeseries.csv", sep=""), sep=";", row.names = F)
    }

    if (save) {
    jpeg(paste(wd,"/output/",main,"_Timeseries.jpg",sep=""), res = 300, width = 8, height = 7, units = 'in')
    par(mfrow=c(1,1),oma=c(1,1,1,1), mgp=c(2, 1,0))
    plot(timeseries[,1],  timeseries[,2], type="b", col="black", pch=16, xlab="year", ylim=c(0,max_index*1.2), ylab=dep_text, main=main.lab) # ylim=c(0,max_index*1.2)
    lines(timeseries[,1], (timeseries[,2]-1.96*timeseries[,3]), type="l",lty=2, col="red" )
    lines(timeseries[,1], (timeseries[,2]+1.96*timeseries[,3]), type="l",lty=2, col="red" )
    legend("topright", c("time series", "CI"), lty=c(1,1), pch=c(16, NA), col=c("black","red"))
    dev.off()
    } else {

        par(mfrow=c(1,1),oma=c(1,1,1,1), mgp=c(2, 1,0))
        plot(timeseries[,1],  timeseries[,2], type="b", col="black", pch=16, xlab="year", ylim=c(0,max_index*1.2), ylab=dep_text, main=main.lab) # ylim=c(0,max_index*1.2)
        lines(timeseries[,1], (timeseries[,2]-1.96*timeseries[,3]), type="l",lty=2, col="red" )
        lines(timeseries[,1], (timeseries[,2]+1.96*timeseries[,3]), type="l",lty=2, col="red" )
        legend("topright", c("time series", "CI"), lty=c(1,1), pch=c(16, NA), col=c("black","red"))
    }













}
