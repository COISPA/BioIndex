#' Sex ratio
#'
#' @author Walter Zupa \email{zupa@fondazionecoispa.org}
#' @param mTATB data frame of the merged TA and TB
#' @param GSA reference GSA for the analysis
#' @param country vector of reference countries for the analysis
#' @param depth_range range of depth strata to perform the analysis (min, max)
#' @param stratas data frame of the reference strata for the study area
#' @param stratification data frame of strata surface area
#' @param wd working directory
#' @param save boolean. If TRUE the plot is saved in the user defined working directory (wd)
#' @param verbose boolean. If TRUE a message is printed
#' @export
sex_ratio <- function(mTATB, GSA, country, depth_range, stratas, stratification, wd=NA, save=TRUE,verbose=FALSE) {

    if (FALSE) {
        GSA=18
        save=TRUE
        merge_TATB <- mTATB
        depth_range <- c(10,800)
        wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\R_BioIndex_3.3_(in update)"
        ta <- read.table(paste(wd,"input/TA/TA GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
        tb <- read.table(paste(wd,"input/TB/TB GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
        tc <- read.table(paste(wd,"input/TC/TC GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
        country="all"
        species <- "MERLMER"
        m <- merge_TATBTC(ta, tb, tc, species="MERLMER", country="all", wd=wd, verbose=TRUE)
        mTATB <- m[[1]]
        sex_ratio(mTATB, GSA, country="all", depth_range, stratification, wd, save)


        stratas <- strata
        stratification <- stratification_tab
    }

    if (is.na(wd) & save) {
        save =FALSE
        if (verbose){
            message("Missing working directory. Results are not saved in the local folder.")
        }
    }

    merge_TATB <- mTATB

    GENERE <- as.character(unique(merge_TATB$GENUS)[unique(merge_TATB$GENUS) != -1])
    SPECIE <- as.character(unique(merge_TATB$SPECIES)[unique(merge_TATB$SPECIES) != -1])
    sspp <- paste(GENERE,SPECIE, sep="")
    GSA <- unique(merge_TATB$GSA)
    species <- sspp

    countries <- unique(mTATB[!is.na(mTATB$COUNTRY),"COUNTRY"])
    l_country <- length(unique(mTATB$COUNTRY))

    if (l_country==1){
        check_country ="Y"
        country_analysis <- countries
    } else {
        if (any(country %in% "all")) {
            country_analysis <<- countries
        } else {
            country_analysis <<- country
        }
    }

    ddd <- merge_TATB
    ddd$NB_FM <- ddd$NB_OF_FEMALES + ddd$NB_OF_MALES
    year_range <- data.frame(year = sort(unique(merge_TATB$YEAR)))

    strata_scheme <- stratas
    strata_scheme <- strata_scheme[strata_scheme$GSA == GSA & strata_scheme$COUNTRY %in% as.character(unique(country_analysis))[1], ]

    depth <- data.frame(matrix(NA,ncol=4, nrow=length(strata_scheme$CODE)))
    colnames(depth) <- c("strata", "min", "max", "bul")
    depth$strata <- strata_scheme$CODE # c(1,2,3,4,5)
    depth$min <-strata_scheme$MIN_DEPTH  # c(10,50,100,200,500)
    depth$max <-strata_scheme$MAX_DEPTH  # c(50,100,200,500,800)
    if (depth_range[2] != 800) {depth_range[2] <- depth_range[2]}
    data <-  ddd[ddd$MEAN_DEPTH>depth_range[1] & ddd$MEAN_DEPTH<=depth_range[2],]
    strata_range <- seq(which(depth[,"min"]==depth_range[1]),which(depth[,"max"]==depth_range[2]),1)
    depth[depth$strata %in% strata_range, "bul"] <- T
    depth[!(depth$strata %in% strata_range), "bul"] <- F

    positive_hauls<- 0

    res_table_calaF <- data.frame(year=sort(unique(merge_TATB$YEAR)))
    res_table_calaF$index_1 <- NA
    res_table_calaF$index_2 <- NA
    res_table_calaF$index_3 <- NA
    res_table_calaF$index_4 <- NA
    res_table_calaF$index_5 <- NA
    res_table_calaF$index_6 <- NA

    res_table_calaFM <- data.frame(year=sort(unique(merge_TATB$YEAR)))
    res_table_calaFM$index_1 <- NA
    res_table_calaFM$index_2 <- NA
    res_table_calaFM$index_3 <- NA
    res_table_calaFM$index_4 <- NA
    res_table_calaFM$index_5 <- NA
    res_table_calaFM$index_6 <- NA

    # FM <- data.frame(year=sort(unique(merge_TATB$YEAR)))
    # FM$index_1 <- NA
    # FM$index_2 <- NA
    # FM$index_3 <- NA
    # FM$index_4 <- NA
    # FM$index_5 <- NA

    SR <- data.frame(year=sort(unique(merge_TATB$YEAR)))
    SR$Indices_F <- NA
    SR$Indices_FM <- NA
    SR$sr <- NA
    SR$variance <- NA

    # stratification <- read.table(paste(wd,"/scripts/utilities/stratification scheme.csv", sep=""), sep=";", header=T)
    stratification$SURF <- as.numeric(stratification$SURF)

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

    res_table_cala2 <- data.frame(year=sort(unique(merge_TATB$YEAR)))
    res_table_cala2$index_1 <- NA
    res_table_cala2$index_2 <- NA
    res_table_cala2$index_3 <- NA
    res_table_cala2$index_4 <- NA
    res_table_cala2$index_5 <- NA
    res_table_cala2$index_6 <- NA

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
    for (i in 1:length(res_table_cala2[,1])) {

        data2 <- ddd[ddd$YEAR == res_table_cala2[i,1] , ]

        if (analysis_stratum1 == T) {s1 <- data2[data2$MEAN_DEPTH >= depth[1,2] & data2$MEAN_DEPTH < depth[1,3], ]}
        if (analysis_stratum2 == T) {s2 <- data2[data2$MEAN_DEPTH >= depth[2,2] & data2$MEAN_DEPTH < depth[2,3], ]}
        if (analysis_stratum3 == T) {s3 <- data2[data2$MEAN_DEPTH >= depth[3,2] & data2$MEAN_DEPTH < depth[3,3], ]}
        if (analysis_stratum4 == T) {s4 <- data2[data2$MEAN_DEPTH >= depth[4,2] & data2$MEAN_DEPTH < depth[4,3], ]}
        if (analysis_stratum5 == T) {s5 <- data2[data2$MEAN_DEPTH >= depth[5,2] & data2$MEAN_DEPTH < depth[5,3], ]}
        if (analysis_stratum6 == T) {s6 <- data2[data2$MEAN_DEPTH >= depth[6,2] & data2$MEAN_DEPTH < depth[6,3], ]}

        if (analysis_stratum1 == T){
            s_nF <- sum(s1[!is.na(s1$NB_OF_FEMALES),"NB_OF_FEMALES"])
            s_nFM <- sum(s1[!is.na(s1$NB_FM),"NB_FM"])
            s_a <- sum(s1$SWEPT_AREA)
            res_table_calaF[i,2] <- s_nF/s_a
            res_table_calaFM[i,2] <- s_nFM/s_a
        } else {area_s1 = 0}

        if (analysis_stratum2 == T){
            s_nF <- sum(s2[!is.na(s2$NB_OF_FEMALES),"NB_OF_FEMALES"])
            s_nFM <- sum(s2[!is.na(s2$NB_FM),"NB_FM"])
            s_a <- sum(s2$SWEPT_AREA)
            res_table_calaF[i,3] <- s_nF/s_a
            res_table_calaFM[i,3] <- s_nFM/s_a
        } else {area_s2 = 0}

        if (analysis_stratum3 == T){
            s_nF <- sum(s3[!is.na(s3$NB_OF_FEMALES),"NB_OF_FEMALES"])
            s_nFM <- sum(s3[!is.na(s3$NB_FM),"NB_FM"])
            s_a <- sum(s3$SWEPT_AREA)
            res_table_calaF[i,4] <- s_nF/s_a
            res_table_calaFM[i,4] <- s_nFM/s_a
        } else {area_s3 = 0}

        if (analysis_stratum4 == T){
            s_nF <- sum(s4[!is.na(s4$NB_OF_FEMALES),"NB_OF_FEMALES"])
            s_nFM <- sum(s4[!is.na(s4$NB_FM),"NB_FM"])
            s_a <- sum(s4$SWEPT_AREA)
            res_table_calaF[i,5] <- s_nF/s_a
            res_table_calaFM[i,5] <- s_nFM/s_a
        } else {area_s4 = 0}

        if (analysis_stratum5 == T){
            #index computation
            s_nF <- sum(s5[!is.na(s5$NB_OF_FEMALES),"NB_OF_FEMALES"])
            s_nFM <- sum(s5[!is.na(s5$NB_FM),"NB_FM"])
            s_a <- sum(s5$SWEPT_AREA)
            res_table_calaF[i,6] <- s_nF/s_a
            res_table_calaFM[i,6] <- s_nFM/s_a
        } else {area_s5 = 0}

        if (analysis_stratum6 == T){
            #index computation
            s_nF <- sum(s6[!is.na(s6$NB_OF_FEMALES),"NB_OF_FEMALES"])
            s_nFM <- sum(s6[!is.na(s6$NB_FM),"NB_FM"])
            s_a <- sum(s6$SWEPT_AREA)
            res_table_calaF[i,7] <- s_nF/s_a
            res_table_calaFM[i,7] <- s_nFM/s_a
        } else {area_s6 = 0}

        positive_hauls[i] <- length(merge_TATB[merge_TATB$TOTAL_NUMBER_IN_THE_HAUL>0 & merge_TATB$YEAR == res_table_cala2[i,1],"TOTAL_NUMBER_IN_THE_HAUL"])

        sum_res_est_F <- c(res_table_calaF[i,2]*peso_s1,res_table_calaF[i,3]*peso_s2,res_table_calaF[i,4]*peso_s3,res_table_calaF[i,5]*peso_s4,res_table_calaF[i,6]*peso_s5,res_table_calaF[i,7]*peso_s6)
        sum_res_est_FM <- c(res_table_calaFM[i,2]*peso_s1,res_table_calaFM[i,3]*peso_s2,res_table_calaFM[i,4]*peso_s3,
                            res_table_calaFM[i,5]*peso_s4,res_table_calaFM[i,6]*peso_s5,res_table_calaFM[i,7]*peso_s6)

        res_table_calaF[i, 8]<- sum(sum_res_est_F[!is.na(sum_res_est_F)])# /sum(area_s1,area_s2,area_s3,area_s4,area_s5)
        res_table_calaFM[i, 8]<- sum(sum_res_est_FM[!is.na(sum_res_est_FM)])# /sum(area_s1,area_s2,area_s3,area_s4,area_s5)
        colnames(res_table_calaF) <- c("year", "stratum 1","stratum 2", "stratum 3", "stratum 4", "stratum 5", "stratum 6", "Indices_F")
        colnames(res_table_calaFM) <- c("year", "stratum 1","stratum 2", "stratum 3", "stratum 4", "stratum 5", "stratum 6", "Indices_FM")
        SR$year[i] <- res_table_calaF$year[i]
        SR$Indices_F[i] <- res_table_calaF$Indices_F[i]
        SR$Indices_FM[i] <- res_table_calaFM$Indices_FM[i]
        SR$sr[i] <- SR$Indices_F[i] / SR$Indices_FM[i]
        SR$variance[i] <- sqrt(SR$sr[i]*(1-SR$sr[i]))/res_table_calaFM[i,8]
    }

    timeseries <- SR


    # plot sex-ratio
    main <- paste(sspp,"_GSA",GSA,"_(Sex ratio)-Random_Stratified_Sampling_",depth_range[1],"-",depth_range[2], " m", sep="")
    main.lab <- paste(sspp," GSA",GSA," (Sex ratio)-RSS ",depth_range[1],"-",depth_range[2], " m", sep="")
    sd_sr <- sqrt(timeseries[!is.na(timeseries$variance),"variance"])
    max_index <- 1 + (max(sd_sr[!is.na(sd_sr)])*1.2)

    if (save){
        write.table(timeseries, paste(wd,"/output/",main,"_Timeseries.csv", sep=""), sep=";", row.names = F)

        jpeg(paste(wd,"/output/",main,"_Timeseries.jpg",sep=""), res = 300, width = 8, height = 7, units = 'in')
        par(mfrow=c(1,1),oma=c(1,1,1,1), mgp=c(2, 1,0))
        plot(timeseries[,1],  timeseries[,4], type="b", col="black", pch=16, xlab="year", ylim=c(0,max_index*1.2), ylab="sex ratio (FF/(FF+MM))", main=main.lab) # ylim=c(0,max_index*1.2)
        lines(timeseries[,1], (timeseries[,4]-1.96*sqrt(timeseries[,5])), type="l",lty=2, col="red" )
        lines(timeseries[,1], (timeseries[,4]+1.96*sqrt(timeseries[,5])), type="l",lty=2, col="red" )
        legend("topright", c("Sex ratio"), lty=c(1,2), pch=c(16), col=c("black"))
        dev.off()
    } else {
        par(mfrow=c(1,1),oma=c(1,1,1,1), mgp=c(2, 1,0))
        plot(timeseries[,1],  timeseries[,4], type="b", col="black", pch=16, xlab="year", ylim=c(0,max_index*1.2), ylab="sex ratio (FF/(FF+MM))", main=main.lab) # ylim=c(0,max_index*1.2)
        lines(timeseries[,1], (timeseries[,4]-1.96*sqrt(timeseries[,5])), type="l",lty=2, col="red" )
        lines(timeseries[,1], (timeseries[,4]+1.96*sqrt(timeseries[,5])), type="l",lty=2, col="red" )
        legend("topright", c("Sex ratio"), lty=c(1,2), pch=c(16), col=c("black"))
    }

return(timeseries)
}
