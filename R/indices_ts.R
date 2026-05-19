#' Estimation of abundance and biomass indices
#' @description
#' Computes time series of stratified abundance and biomass indices,
#' serving as the primary tool for monitoring temporal trends in stock status.
#'
#' @author Walter Zupa \email{zupa@fondazionecoispa.org}
#' @param mTATB data frame
#' @param GSA reference GSA for the analysis
#' @param country reference countries in the GSA for the analysis
#' @param depth_range range of depth strata to perform the analysis (min, max)
#' @param strata_scheme data frame of the stratification scheme
#' @param stratification data frame of strata surface area
#' @param wd working directory
#' @param save boolean. If TRUE the plot is saved in the user defined working directory (wd)
#' @importFrom gridExtra grid.arrange
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_density ggtitle geom_boxplot xlab ylab ylim xlim ggsave theme_bw geom_line geom_point
#' @importFrom stringr str_split
#' @import grDevices
#' @import graphics
#' @importFrom utils write.table
#' @importFrom grDevices jpeg dev.off
#' @importFrom stats sd
#' @return A \code{list} of \code{ggplot} objects representing the time series of calculated indices.
#' @examples
#' data(TA)
#' data(TB)
#' # Create a merged dataset for GSA 10
#' mTATB <- merge_TATB(TA[TA$AREA==10,], TB[TB$AREA==10,], "MERLMER")
#' # Run indices_ts
#' indices_ts(mTATB, GSA=10, country="all", depth_range=c(10,800),
#'            strata_scheme=BioIndex::strata_scheme,
#'            stratification=BioIndex::stratification, wd=tempdir(), save=FALSE)
#' @export

indices_ts <- function(mTATB, GSA, country="all", depth_range, strata_scheme, stratification, wd=NA, save=TRUE) {
    if (is.na(wd)) { wd <- tempdir() }
    oldpar <- par(no.readonly = TRUE)
    oldpar$new <- NULL
    on.exit(suppressWarnings(par(oldpar)))

    if (FALSE) {

        GSA=10
        save=TRUE
        merge_TATB <- mTATB
        depth_range <- c(10,800)

        wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\R_BioIndex_3.3_(in update)"
        # strata_scheme <- read.table(paste(wd, "/scripts/utilities/strata.csv", sep=""), sep=";",header = T)
        # save(strata_scheme,file="D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\- Pacchetto R BioIndex -\\BioIndex\\data\\strata_scheme.rda",compress="xz",compression_level=9)

        # stratification <- read.table(paste(wd,"/scripts/utilities/stratification scheme.csv", sep=""), sep=";", header=T)
        # save(stratification,file="D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\- Pacchetto R BioIndex -\\BioIndex\\data\\stratification.rda",compress="xz",compression_level=9)

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
        species <- "MERLMER"

        m <- merge_TATBTC(ta, tb, tc, species="MERLMER", country="all", wd=wd, verbose=TRUE)
        mTATB <- m[[1]]
        indices_ts(mTATB, GSA=10, country="all", depth_range=c(10,800), strata_scheme=BioIndex::strata_scheme, stratification=BioIndex::stratification,wd, save=TRUE)

        # indices_ts(mTATB, GSA=10, country="all", depth_range=c(10,800), strata_scheme,wd, stratification, save=TRUE)
    }

  if (is.na(wd) & save) {
    save=TRUE
    if (verbose){
      message("Missing working directory. Results are not saved in the local folder.")
    }
  }

    MEAN_DEPTH <- strata <- value <- year <- NULL

    merge_TATB <- mTATB
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

    merge_TATB <- merge_TATB[merge_TATB$COUNTRY %in% country_analysis, ]

    strata_scheme <- strata_scheme[strata_scheme$GSA == GSA & strata_scheme$COUNTRY %in% as.character(unique(country_analysis))[1], ]

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
    # plotting species' depth range information
    g1 <- ggplot(merge_TATB[merge_TATB[,varcol] >0,], aes(x=MEAN_DEPTH)) + geom_density()+xlab("depth (m)")+ggtitle(paste(sspp, " - GSA",GSA,sep="" ))
    if (save){
        if(save) ggsave(paste(wd,"/output/depth.distribution_",sspp,"_GSA",GSA,".jpg", sep=""),width=5, height=5, plot = g1)
        }
    g2 <- ggplot(merge_TATB[merge_TATB[,varcol]>0,], aes(x=species, y=MEAN_DEPTH)) + geom_boxplot() + xlab("") + ylim(0, max(merge_TATB$MEAN_DEPTH))+ylab("depth (m)")+ggtitle(paste(sspp, " - GSA",GSA,sep="" ))
    if (save){
        if(save) ggsave(paste(wd,"/output/depth.distribution_(boxplot)",sspp,"_GSA",GSA,".jpg", sep=""),width=5, height=5, plot = g2)
    }
    grid.arrange(g1, g2, ncol=2)
    #------------------------------------------

    #------------------------------------------
    # selection of depth range
    data <- merge_TATB
    depth <- data.frame(matrix(NA,ncol=4, nrow=length(strata_scheme$CODE)))
    colnames(depth) <- c("strata", "min", "max", "bul")
    depth$strata <- strata_scheme$CODE # c(1,2,3,4,5)
    depth$min <-strata_scheme$MIN_DEPTH  # c(10,50,100,200,500)
    depth$max <-strata_scheme$MAX_DEPTH  # c(50,100,200,500,800)
    # depth_range <- dlgInput("Depth range for the analysis: ", default="5,35",Sys.info()[""])$res
    # depth_range <- data.frame(range = strsplit(depth_range, ",")); depth_range <- as.numeric(as.character(depth_range[,1]))
    if (depth_range[2] != 800) {depth_range[2] <- depth_range[2]}
    data <-  data[data$MEAN_DEPTH>depth_range[1] & data$MEAN_DEPTH<=depth_range[2],]
    idx_from <- which(depth[,"min"] == depth_range[1])
    idx_to   <- which(depth[,"max"] == depth_range[2])

    if (length(idx_from) != 1 || length(idx_to) != 1) {
        stop(paste0("Error: depth range [", depth_range[1], ", ", depth_range[2],
                    "] does not match the strata boundaries defined for this GSA. ",
                    "Available min depths: ", paste(sort(unique(depth$min)), collapse=", "),
                    "; Available max depths: ", paste(sort(unique(depth$max)), collapse=", ")))
    }
    strata_range <- seq(idx_from, idx_to, 1)
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
    years = sort(unique(merge_TATB$YEAR))
    res_table_cala2 <- data.frame(year=sort(unique(merge_TATB$YEAR)))
    res_table_cala2$index_1 <- NA
    res_table_cala2$index_2 <- NA
    res_table_cala2$index_3 <- NA
    res_table_cala2$index_4 <- NA
    res_table_cala2$index_5 <- NA
    res_table_cala2$index_6 <- NA

    sd_strata <- data.frame(year=sort(unique(merge_TATB$YEAR)))
    sd_strata$sd_1 <- NA
    sd_strata$sd_2 <- NA
    sd_strata$sd_3 <- NA
    sd_strata$sd_4 <- NA
    sd_strata$sd_5 <- NA
    sd_strata$sd_6 <- NA

    se_table2  <-  data.frame(year=sort(unique(merge_TATB$YEAR)))
    se_table2$sum_1 <- NA
    se_table2$sum_2 <- NA
    se_table2$sum_3 <- NA
    se_table2$sum_4 <- NA
    se_table2$sum_5 <- NA
    se_table2$sum_6 <- NA
    se_table2$sd <- NA

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

    i=15
    for(i in 1:length(res_table_cala2[,1])) {

        data2 <- ddd[ddd$YEAR == res_table_cala2[i,1] , ]

        if (analysis_stratum1 == T) {s1 <- data2[floor(data2$MEAN_DEPTH) >= depth[1,2] & floor(data2$MEAN_DEPTH) <= depth[1,3], ]}
        if (analysis_stratum2 == T) {s2 <- data2[floor(data2$MEAN_DEPTH) > depth[2,2] & floor(data2$MEAN_DEPTH) <= depth[2,3], ]}
        if (analysis_stratum3 == T) {s3 <- data2[floor(data2$MEAN_DEPTH) > depth[3,2] & floor(data2$MEAN_DEPTH) <= depth[3,3], ]}
        if (analysis_stratum4 == T) {s4 <- data2[floor(data2$MEAN_DEPTH) > depth[4,2] & floor(data2$MEAN_DEPTH) <= depth[4,3], ]}
        if (analysis_stratum5 == T) {s5 <- data2[floor(data2$MEAN_DEPTH) > depth[5,2] & floor(data2$MEAN_DEPTH) <= depth[5,3], ]}
        if (analysis_stratum6 == T) {s6 <- data2[floor(data2$MEAN_DEPTH) > depth[6,2] & floor(data2$MEAN_DEPTH) <= depth[6,3], ]}

        if (analysis_stratum1 == T){
            s_n <- sum(s1[!is.na(s1$TOTAL_NUMBER_IN_THE_HAUL),"TOTAL_NUMBER_IN_THE_HAUL"])
            s_a <- sum(s1$SWEPT_AREA)
            res_table_cala2[i,2] <- s_n/s_a

            mean_ind <- mean(s1[!is.na(s1$N_km2), "N_km2" ])
            sq <- (s1[!is.na(s1$N_km2),"N_km2"] - mean_ind)^2
            sqA <- sq*s1$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s1[,1])-1)
            f1 <- s_a/area_s1
            se_table2[i,2] <- (((peso_s1)^2 * sum_sqA)/ s_a)* (1-f1)
            sd_strata[i,2] <- sqrt( ((sum_sqA)/ s_a)* (1-f1))
        } else {area_s1 = 0}

        if (analysis_stratum2 == T){
            s_n <- sum(s2[!is.na(s2$TOTAL_NUMBER_IN_THE_HAUL),"TOTAL_NUMBER_IN_THE_HAUL"])
            s_a <- sum(s2$SWEPT_AREA)
            res_table_cala2[i,3] <- s_n/s_a

            mean_ind <- mean(s2[!is.na(s2$N_km2), "N_km2" ])
            sq <- (s2[!is.na(s2$N_km2),"N_km2"] - mean_ind)^2
            sqA <- sq*s2$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s2[,1])-1)
            f2 <- s_a/area_s2
            se_table2[i,3] <- (((peso_s2)^2 * sum_sqA)/ s_a)* (1-f2)
            sd_strata[i,3] <- sqrt( ((sum_sqA)/ s_a)* (1-f2))
        } else {area_s2 = 0}

        if (analysis_stratum3 == T){
            s_n <- sum(s3[!is.na(s3$TOTAL_NUMBER_IN_THE_HAUL),"TOTAL_NUMBER_IN_THE_HAUL"])
            s_a <- sum(s3$SWEPT_AREA)
            res_table_cala2[i,4] <- s_n/s_a

            mean_ind <- mean(s3[!is.na(s3$N_km2), "N_km2" ])
            sq <- (s3[!is.na(s3$N_km2),"N_km2"] - mean_ind)^2
            sqA <- sq*s3$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s3[,1])-1)
            f3 <- s_a/area_s3
            se_table2[i,4] <- (((peso_s3)^2 * sum_sqA)/ s_a)* (1-f3)
            sd_strata[i,4] <- sqrt( ((sum_sqA)/ s_a)* (1-f3))
        } else {area_s3 = 0}

        if (analysis_stratum4 == T){
            s_n <- sum(s4[!is.na(s4$TOTAL_NUMBER_IN_THE_HAUL),"TOTAL_NUMBER_IN_THE_HAUL"])
            s_a <- sum(s4$SWEPT_AREA)
            res_table_cala2[i,5] <- s_n/s_a

            mean_ind <- mean(s4[!is.na(s4$N_km2), "N_km2" ])
            sq <- (s4[!is.na(s4$N_km2),"N_km2"] - mean_ind)^2
            sqA <- sq*s4$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s4[,1])-1)
            f4 <- s_a/area_s4
            se_table2[i,5] <- (((peso_s4)^2 * sum_sqA)/ s_a)* (1-f4)
            sd_strata[i,5] <- sqrt( ((sum_sqA)/ s_a)* (1-f4))
        } else {area_s4 = 0}

        if (analysis_stratum5 == T){
            #index computation
            s_n <- sum(s5[!is.na(s5$TOTAL_NUMBER_IN_THE_HAUL),"TOTAL_NUMBER_IN_THE_HAUL"])
            s_a <- sum(s5$SWEPT_AREA)
            res_table_cala2[i,6] <- s_n/s_a

            #se computation
            mean_ind <- mean(s5[!is.na(s5$N_km2), "N_km2" ])
            sq <- (s5[!is.na(s5$N_km2),"N_km2"] - mean_ind)^2
            sqA <- sq*s5$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s5[,1])-1)
            f5 <- s_a/area_s5
            se_table2[i,6] <- (((peso_s5)^2 * sum_sqA)/ s_a)* (1-f5)
            sd_strata[i,6] <- sqrt( ((sum_sqA)/ s_a)* (1-f5))
        } else {area_s5 = 0}

        if (analysis_stratum6 == T){
            #index computation
            s_n <- sum(s6[!is.na(s6$TOTAL_NUMBER_IN_THE_HAUL),"TOTAL_NUMBER_IN_THE_HAUL"])
            s_a <- sum(s6$SWEPT_AREA)
            res_table_cala2[i,7] <- s_n/s_a

            #se computation
            mean_ind <- mean(s6[!is.na(s6$N_km2), "N_km2" ])
            sq <- (s6[!is.na(s6$N_km2),"N_km2"] - mean_ind)^2
            sqA <- sq*s6$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s6[,1])-1)
            f6 <- s_a/area_s6
            se_table2[i,7] <- (((peso_s6)^2 * sum_sqA)/ s_a)* (1-f6)
            sd_strata[i,7] <- sqrt( ((sum_sqA)/ s_a)* (1-f6))
        } else {area_s6 = 0}

        # source(paste(wd, "/scripts/Index_estimation_cala.R", sep=""), encoding = 'UTF-8')
        positive_hauls[i] <- length(merge_TATB[merge_TATB$TOTAL_NUMBER_IN_THE_HAUL>0 & merge_TATB$YEAR == res_table_cala2[i,1],"TOTAL_NUMBER_IN_THE_HAUL"])
        total_hauls[i] <- length(merge_TATB[merge_TATB$YEAR == res_table_cala2[i,1],"HAUL_NUMBER"])

        sum_res_est <- c(res_table_cala2[i,2]*peso_s1,res_table_cala2[i,3]*peso_s2,res_table_cala2[i,4]*peso_s3,res_table_cala2[i,5]*peso_s4,res_table_cala2[i,6]*peso_s5,res_table_cala2[i,7]*peso_s6)
        res_table_cala2[i, 8]<- sum(sum_res_est[!is.na(sum_res_est)])# /sum(area_s1,area_s2,area_s3,area_s4,area_s5)
        colnames(res_table_cala2) <- c("year", "stratum 1","stratum 2", "stratum 3", "stratum 4", "stratum 5", "stratum 6", "Indices")
        se_table2[i, "sd"] <- sqrt(rowSums(se_table2[i, 2:7], na.rm = T))
    }

    timeseries <- data.frame(year = res_table_cala2[,1], abundance = res_table_cala2[,8], sd= se_table2[, "sd"])
    timeseries$CV <- timeseries[,3] / timeseries[,2]
    timeseries$invCV <- 1/timeseries$CV
    timeseries$positive_hauls_perc <- positive_hauls / total_hauls * 100
    colnames(res_table_cala2) <- c("year", "stratum 1","stratum 2", "stratum 3", "stratum 4", "stratum 5", "stratum 6", "Indices")

    {
        cols.strata <- (depth[!is.na(depth$strata), 'strata'] + 1)
        n.col.strata <- length(cols.strata)
        cv_strata <- sd_strata
        col=5
        for (col in cols.strata) {
            row=1
            for (row in 1:length(year_range[,1])) {
                if (!is.na(sd_strata[row,col])) {
                    if (res_table_cala2[row,col]!=0){
                        cv_strata[row,col] <- cv_strata[row,col]/res_table_cala2[row,col]
                    } else {
                        cv_strata[row,col] <- 0
                    }
                }
            }
        }

        TAB_SD_CV <- cbind(sd_strata[, c(1,cols.strata)],cv_strata[,c(cols.strata)])
        colnames(TAB_SD_CV) <- c("year", paste("SD", depth[!is.na(depth$strata), 'strata'],sep="_"), paste("CV", depth[!is.na(depth$strata), 'strata'],sep="_"))
        tab.strata.abu <- data.frame(year=res_table_cala2[,"year"], res_table_cala2[, cols.strata], TAB_SD_CV[,-1])
        colnames(tab.strata.abu) <- c("year", paste("Abundance", depth[!is.na(depth$strata), 'strata'],sep="_"), paste("SD", depth[!is.na(depth$strata), 'strata'],sep="_"), paste("CV", depth[!is.na(depth$strata), 'strata'],sep="_"))

        title.abu.haul <- paste(sspp,"_GSA",GSA,"_(",dependent,"_by_stratum)-RSS_",depth_range[1],"-",depth_range[2], " m", sep="")

        tab.strata.abu
        if (save){
        write.table(tab.strata.abu, paste(wd,"/output/",title.abu.haul,".csv", sep=""), sep=";", row.names = F)
        }

        sel.strata1 <- c(1,depth[(depth$bul), "strata"] + 1)
        sel.strata2 <- c(1, depth[(depth$bul), "strata"] + 1 + length(cols.strata))
        tab.ind <- tab.strata.abu[, sel.strata1]
        tab.sd <- sd_strata[, sel.strata1]
        df.plot1 <- melt(tab.ind, id.vars = 1)
        df.plot2 <- melt(tab.sd, id.vars = 1)
        var1 <-  as.data.frame(do.call(rbind, str_split(as.character(df.plot1$variable),"_"))); colnames(var1) <- c("var", "strata")
        var2 <-  as.data.frame(do.call(rbind, str_split(as.character(df.plot2$variable),"_"))); colnames(var2) <- c("var", "strata")
        df.plot1 <- cbind(df.plot1,var1)
        df.plot2 <- cbind(df.plot2,var2)
        df.plot1 <- data.frame(df.plot1, sd=df.plot2$value)
        main <- paste(sspp,"_GSA",GSA,"_(",dependent,"_by_stratum)-RSS_",depth_range[1],"-",depth_range[2], " m", sep="")
        main.lab <- paste(sspp," GSA",GSA," (",dependent," by stratum)-RSS ",depth_range[1],"-",depth_range[2], " m", sep="")

        if (length(years)>1){
            p1 <- ggplot(data=df.plot1, aes(x=year, y=value, group=strata, colour=strata)) +
                geom_line() +
                geom_point() +
                xlab("year") +
                ylab(dep_text) +
                ggtitle(main.lab) +
                theme_bw()
        # print(p1) # rimosso per policy CRAN
            if(save){
                if(save) ggsave(paste(wd,"/output/",main,"_Timeseries.jpg",sep=""), dpi=300 , width=6, height=5, plot = p1) #
                }
        } else {
            p1 <- ggplot(data=df.plot1, aes(x=year, y=value, group=strata, colour=strata)) +
                geom_point() +
                xlab("year") +
                ylab(dep_text) +
                ggtitle(main.lab) +
                theme_bw()
        # print(p1) # rimosso per policy CRAN
            if(save){
                if(save) ggsave(paste(wd,"/output/",main,"_Timeseries.jpg",sep=""), dpi=300 , width=6, height=5, plot = p1)
            }
        }

    }


    # plot timeseries of mean indices
    main <- paste(sspp,"_GSA",GSA,"_(",dependent,")-Random_Stratified_Sampling_",depth_range[1],"-",depth_range[2], " m", sep="")
    main.lab <- paste(sspp," GSA",GSA," (",dependent,")-RSS ",depth_range[1],"-",depth_range[2], " m", sep="")
    max_index <- max(timeseries[,2]) + (max(timeseries[!is.na(timeseries$sd),"sd"])*1.2)

    if(save){
        write.table(timeseries, paste(wd,"/output/",main,"_Timeseries.csv", sep=""), sep=";", row.names = F)
}
timeseries_abundance <- timeseries
    if (save){
    jpeg(filename =paste(wd,"/output/",main,"_Timeseries.jpg",sep=""), res = 300, width = 8, height = 7, units = 'in')
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


    # plot timeseries of inverse of CV
    main <- paste(sspp,"_GSA",GSA,"_(inverseCV of abundance)-Random_Stratified_Sampling_",depth_range[1],"-",depth_range[2], " m", sep="")
    main.lab <- paste(sspp," GSA",GSA," (1/CV of mean abundance)-RSS ",depth_range[1],"-",depth_range[2], " m", sep="")
    max_index <- max(timeseries[!is.na(timeseries$invCV) & timeseries$invCV != Inf,"invCV"])

  if (save){
    jpeg(filename =paste(wd,"/output/",main,"_Timeseries.jpg",sep=""), res = 300, width = 8, height = 7, units = 'in')
    par(mfrow=c(1,1),oma=c(1,1,1,1), mgp=c(2, 1,0))
    plot(timeseries[,1],  timeseries[,5], type="b", col="black", pch=16, xlab="year", ylim=c(0,max_index*1.2), ylab="1/CV", main=main.lab) # ylim=c(0,max_index*1.2)
    legend("topright", c("inverse CV"), lty=c(1), pch=c(16), col=c("black"))
    dev.off()
  } else {
      par(mfrow=c(1,1),oma=c(1,1,1,1), mgp=c(2, 1,0))
      plot(timeseries[,1],  timeseries[,5], type="b", col="black", pch=16, xlab="year", ylim=c(0,max_index*1.2), ylab="1/CV", main=main.lab) # ylim=c(0,max_index*1.2)
      legend("topright", c("inverse CV"), lty=c(1), pch=c(16), col=c("black"))
  }

    # plot number of positive hauls per year
    main <- paste(sspp,"_GSA",GSA,"_(positive hauls)-Random_Stratified_Sampling_",depth_range[1],"-",depth_range[2], " m", sep="")
    main.lab <- paste(sspp," GSA",GSA," (positive hauls %)-RSS ",depth_range[1],"-",depth_range[2], " m", sep="")
    max_index <- max(timeseries[!is.na(timeseries$positive_hauls_perc),"positive_hauls_perc"])
    if (save) {
        jpeg(filename =paste(wd,"/output/",main,"_Timeseries.jpg",sep=""), res = 300, width = 8, height = 7, units = 'in')
        par(mfrow=c(1,1),oma=c(1,1,1,1), mgp=c(2, 1,0))
        plot(timeseries[,1],  timeseries[,6], type="b", col="black", pch=16, xlab="year", ylim=c(0,max_index*1.2), ylab="number of hauls (%)", main=main.lab) # ylim=c(0,max_index*1.2)
        legend("topright", c("number of positive hauls %"), lty=c(1), pch=c(16), col=c("black"))
        dev.off()
    } else {
        par(mfrow=c(1,1),oma=c(1,1,1,1), mgp=c(2, 1,0))
        plot(timeseries[,1],  timeseries[,6], type="b", col="black", pch=16, xlab="year", ylim=c(0,max_index*1.2), ylab="number of hauls (%)", main=main.lab) # ylim=c(0,max_index*1.2)
        legend("topright", c("number of positive hauls %"), lty=c(1), pch=c(16), col=c("black"))
    }
    message("\n Estimation of abundance indices completed ")

    index_ts_M(mTATB, GSA, country_analysis, depth_range, strata_scheme, stratification, wd=wd, save)
    index_ts_F(mTATB, GSA, country_analysis, depth_range, strata_scheme, stratification, wd=wd, save)

    #-----------------------
    # BIOMASS INDICES
    #-----------------------

    #------------------------
    # Variable definitions #
    dependent <- "biomass"
    dep_text <-expression(paste("Biomass ", (kg/km^2), sep=" "))
    varcol <- which(colnames(merge_TATB)=="kg_km2")
    col_response <- which(colnames(merge_TATB)=="kg_km2")
    GENERE <- as.character(unique(merge_TATB$GENUS)[unique(merge_TATB$GENUS) != -1])
    SPECIE <- as.character(unique(merge_TATB$SPECIES)[unique(merge_TATB$SPECIES) != -1])
    sspp <- paste(GENERE,SPECIE, sep="")
    GSA <- unique(merge_TATB$GSA)
    species <- sspp
    #------------------------

    ddd <- merge_TATB
    year_range <- data.frame(year = sort(unique(merge_TATB$YEAR)))

    res_table_cala2 <- data.frame(year=sort(unique(merge_TATB$YEAR)))
    res_table_cala2$index_1 <- NA
    res_table_cala2$index_2 <- NA
    res_table_cala2$index_3 <- NA
    res_table_cala2$index_4 <- NA
    res_table_cala2$index_5 <- NA
    res_table_cala2$index_6 <- NA

    sd_strata <- data.frame(year=sort(unique(merge_TATB$YEAR)))
    sd_strata$sd_1 <- NA
    sd_strata$sd_2 <- NA
    sd_strata$sd_3 <- NA
    sd_strata$sd_4 <- NA
    sd_strata$sd_5 <- NA
    sd_strata$sd_6 <- NA

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

        if (analysis_stratum1 == T) {s1 <- data2[data2$MEAN_DEPTH > depth[1,2] & data2$MEAN_DEPTH <= depth[1,3], ]}
        if (analysis_stratum2 == T) {s2 <- data2[data2$MEAN_DEPTH > depth[2,2] & data2$MEAN_DEPTH <= depth[2,3], ]}
        if (analysis_stratum3 == T) {s3 <- data2[data2$MEAN_DEPTH > depth[3,2] & data2$MEAN_DEPTH <= depth[3,3], ]}
        if (analysis_stratum4 == T) {s4 <- data2[data2$MEAN_DEPTH > depth[4,2] & data2$MEAN_DEPTH <= depth[4,3], ]}
        if (analysis_stratum5 == T) {s5 <- data2[data2$MEAN_DEPTH > depth[5,2] & data2$MEAN_DEPTH <= depth[5,3], ]}
        if (analysis_stratum6 == T) {s6 <- data2[data2$MEAN_DEPTH > depth[6,2] & data2$MEAN_DEPTH <= depth[6,3], ]}

        if (analysis_stratum1 == T){
            s_n <- sum(s1[!is.na(s1$TOTAL_WEIGHT_IN_THE_HAUL),"TOTAL_WEIGHT_IN_THE_HAUL"]/1000)
            s_a <- sum(s1$SWEPT_AREA)
            res_table_cala2[i,2] <- s_n/s_a

            mean_ind <- mean(s1[!is.na(s1$kg_km2), "kg_km2" ])
            sq <- (s1[!is.na(s1$kg_km2),"kg_km2"] - mean_ind)^2
            sqA <- sq*s1$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s1[,1])-1)
            f1 <- s_a/area_s1
            se_table2[i,2] <- (((peso_s1)^2 * sum_sqA)/ s_a)* (1-f1)
            sd_strata[i,2] <- sqrt( ((sum_sqA)/ s_a)* (1-f1))

        } else {area_s1 = 0}

        if (analysis_stratum2 == T){
            s_n <- sum(s2[!is.na(s2$TOTAL_WEIGHT_IN_THE_HAUL),"TOTAL_WEIGHT_IN_THE_HAUL"]/1000)
            s_a <- sum(s2$SWEPT_AREA)
            res_table_cala2[i,3] <- s_n/s_a

            mean_ind <- mean(s2[!is.na(s2$kg_km2), "kg_km2" ])
            sq <- (s2[!is.na(s2$kg_km2),"kg_km2"] - mean_ind)^2
            sqA <- sq*s2$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s2[,1])-1)
            f2 <- s_a/area_s2
            se_table2[i,3] <- (((peso_s2)^2 * sum_sqA)/ s_a)* (1-f2)
            sd_strata[i,3] <- sqrt( ((sum_sqA)/ s_a)* (1-f2))
        } else {area_s2 = 0}

        if (analysis_stratum3 == T){
            s_n <- sum(s3[!is.na(s3$TOTAL_WEIGHT_IN_THE_HAUL),"TOTAL_WEIGHT_IN_THE_HAUL"]/1000)
            s_a <- sum(s3$SWEPT_AREA)
            res_table_cala2[i,4] <- s_n/s_a

            mean_ind <- mean(s3[!is.na(s3$kg_km2), "kg_km2" ])
            sq <- (s3[!is.na(s3$kg_km2),"kg_km2"] - mean_ind)^2
            sqA <- sq*s3$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s3[,1])-1)
            f3 <- s_a/area_s3
            se_table2[i,4] <- (((peso_s3)^2 * sum_sqA)/ s_a)* (1-f3)
            sd_strata[i,4] <- sqrt( ((sum_sqA)/ s_a)* (1-f3))
        } else {area_s3 = 0}

        if (analysis_stratum4 == T){
            s_n <- sum(s4[!is.na(s4$TOTAL_WEIGHT_IN_THE_HAUL),"TOTAL_WEIGHT_IN_THE_HAUL"]/1000)
            s_a <- sum(s4$SWEPT_AREA)
            res_table_cala2[i,5] <- s_n/s_a

            mean_ind <- mean(s4[!is.na(s4$kg_km2), "kg_km2" ])
            sq <- (s4[!is.na(s4$kg_km2),"kg_km2"] - mean_ind)^2
            sqA <- sq*s4$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s4[,1])-1)
            f4 <- s_a/area_s4
            se_table2[i,5] <- (((peso_s4)^2 * sum_sqA)/ s_a)* (1-f4)
            sd_strata[i,5] <- sqrt( ((sum_sqA)/ s_a)* (1-f4))
        } else {area_s4 = 0}

        if (analysis_stratum5 == T){
            #index computation
            s_n <- sum(s5[!is.na(s5$TOTAL_WEIGHT_IN_THE_HAUL),"TOTAL_WEIGHT_IN_THE_HAUL"]/1000)
            s_a <- sum(s5$SWEPT_AREA)
            res_table_cala2[i,6] <- s_n/s_a

            #se computation
            mean_ind <- mean(s5[!is.na(s5$kg_km2), "kg_km2" ])
            sq <- (s5[!is.na(s5$kg_km2),"kg_km2"] - mean_ind)^2
            sqA <- sq*s5$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s5[,1])-1)
            f5 <- s_a/area_s5
            se_table2[i,6] <- (((peso_s5)^2 * sum_sqA)/ s_a)* (1-f5)
            sd_strata[i,6] <- sqrt( ((sum_sqA)/ s_a)* (1-f5))
        } else {area_s5 = 0}

        if (analysis_stratum6 == T){
            #index computation
            s_n <- sum(s6[!is.na(s6$TOTAL_WEIGHT_IN_THE_HAUL),"TOTAL_WEIGHT_IN_THE_HAUL"]/1000)
            s_a <- sum(s6$SWEPT_AREA)
            res_table_cala2[i,7] <- s_n/s_a

            #se computation
            mean_ind <- mean(s6[!is.na(s6$kg_km2), "kg_km2" ])
            sq <- (s6[!is.na(s6$kg_km2),"kg_km2"] - mean_ind)^2
            sqA <- sq*s6$SWEPT_AREA
            sum_sqA <- sum(sqA)/(length(s6[,1])-1)
            f6 <- s_a/area_s6
            se_table2[i,7] <- (((peso_s6)^2 * sum_sqA)/ s_a)* (1-f6)
            sd_strata[i,7] <- sqrt( ((sum_sqA)/ s_a)* (1-f6))
        } else {area_s6 = 0}


        sum_res_est <- c(res_table_cala2[i,2]*peso_s1,res_table_cala2[i,3]*peso_s2,res_table_cala2[i,4]*peso_s3,res_table_cala2[i,5]*peso_s4,res_table_cala2[i,6]*peso_s5,res_table_cala2[i,7]*peso_s6)
        res_table_cala2[i, 8]<- sum(sum_res_est[!is.na(sum_res_est)])# /sum(area_s1,area_s2,area_s3,area_s4,area_s5)
        colnames(res_table_cala2) <- c("year", "stratum 1","stratum 2", "stratum 3", "stratum 4", "stratum 5", "stratum 6", "Indices")
        se_table2[i, "sd"] <-sqrt(rowSums(se_table2[i, 2:7], na.rm = T))
    }

    timeseries <- data.frame(year = res_table_cala2[,1], biomass = res_table_cala2[,8], sd= se_table2[, "sd"])
    timeseries$CV <- timeseries[,3] / timeseries[,2]
    timeseries$invCV <- 1/timeseries$CV
    colnames(res_table_cala2) <- c("year", "stratum 1","stratum 2", "stratum 3", "stratum 4", "stratum 5", "stratum 6", "Indices")


    {
        cols.strata <- (depth[!is.na(depth$strata), 'strata'] + 1)
        n.col.strata <- length(cols.strata)
        cv_strata <- sd_strata
        col=5
        for (col in cols.strata) {
            row=1
            for (row in 1:length(year_range[,1])) {
                if (!is.na(sd_strata[row,col])) {
                    if (res_table_cala2[row,col]!=0){
                        cv_strata[row,col] <- cv_strata[row,col]/res_table_cala2[row,col]
                    } else {
                        cv_strata[row,col] <- 0
                    }
                }
            }
        }


        TAB_SD_CV <- cbind(sd_strata[, c(1,cols.strata)],cv_strata[,c(cols.strata)])
        colnames(TAB_SD_CV) <- c("year", paste("SD", depth[!is.na(depth$strata), 'strata'],sep="_"), paste("CV", depth[!is.na(depth$strata), 'strata'],sep="_"))
        tab.strata.bio <- data.frame(year=res_table_cala2[,"year"], res_table_cala2[, cols.strata], TAB_SD_CV[,-1])
        colnames(tab.strata.bio) <- c("year", paste("Biomass", depth[!is.na(depth$strata), 'strata'],sep="_"), paste("SD", depth[!is.na(depth$strata), 'strata'],sep="_"), paste("CV", depth[!is.na(depth$strata), 'strata'],sep="_"))

        title.bio.haul <- paste(sspp,"_GSA",GSA,"_(",dependent,"_by_stratum)-RSS_",depth_range[1],"-",depth_range[2], " m", sep="")
        if (save){
            write.table(tab.strata.bio, paste(wd,"/output/",title.bio.haul,".csv", sep=""), sep=";", row.names = F)
        }
        sel.strata1 <- c(1,depth[(depth$bul), "strata"] + 1)
        sel.strata2 <- c(1, depth[(depth$bul), "strata"] + 1 + length(cols.strata))
        tab.ind <- tab.strata.bio[, sel.strata1]
        tab.sd <- sd_strata[, sel.strata1]
        df.plot1 <- melt(tab.ind, id.vars = 1)
        df.plot2 <- melt(tab.sd, id.vars = 1)
        var1 <-  as.data.frame(do.call(rbind, str_split(as.character(df.plot1$variable),"_"))); colnames(var1) <- c("var", "strata")
        var2 <-  as.data.frame(do.call(rbind, str_split(as.character(df.plot2$variable),"_"))); colnames(var2) <- c("var", "strata")
        df.plot1 <- cbind(df.plot1,var1)
        df.plot2 <- cbind(df.plot2,var2)
        df.plot1 <- data.frame(df.plot1, sd=df.plot2$value)
        main <- paste(sspp,"_GSA",GSA,"_(",dependent,"_by_stratum)-RSS_",depth_range[1],"-",depth_range[2], " m", sep="")
        main.lab <- paste(sspp," GSA",GSA," (",dependent," by stratum)-RSS ",depth_range[1],"-",depth_range[2], " m", sep="")

        if (length(years)>1){
            p2 <- ggplot(data=df.plot1, aes(x=year, y=value, group=strata, colour=strata)) +
                geom_line() +
                geom_point() +
                xlab("year") +
                ylab(dep_text) +
                ggtitle(main.lab) +
                theme_bw()
        # print(p2) # rimosso per policy CRAN
            if (save) {ggsave(paste(wd,"/output/",main,"_Timeseries.jpg",sep=""), dpi=300 , width=6, height=5, plot = p2)} #
        } else {
           p2 <-  ggplot(data=df.plot1, aes(x=year, y=value, group=strata, colour=strata)) +
                geom_point() +
                xlab("year") +
                ylab(dep_text) +
                ggtitle(main.lab) +
                theme_bw()
        # print(p2) # rimosso per policy CRAN
           if (save) {ggsave(paste(wd,"/output/",main,"_Timeseries.jpg",sep=""), dpi=300 , width=6, height=5, plot = p2)}
        }
    }
    # is.na(timeseries) <- do.call(cbind,lapply(timeseries, is.infinite))
    # plot timeseries of mean indices
    main <- paste(sspp,"_GSA",GSA,"_(",dependent,")-Random_Stratified_Sampling_",depth_range[1],"-",depth_range[2], " m", sep="")
    main.lab <- paste(sspp," GSA",GSA," (",dependent,")-RSS ",depth_range[1],"-",depth_range[2], " m", sep="")
    max_index <- max(timeseries[,2])  + (max(timeseries[!is.na(timeseries$sd),"sd"])*1.2)
    if (save) {write.table(timeseries, paste(wd,"/output/",main,"_Timeseries.csv", sep=""), sep=";", row.names = F)}
    timeseries_biomass <- timeseries
if (save){
    jpeg(filename =paste(wd,"/output/",main,"_Timeseries.jpg",sep=""), res = 300, width = 8, height = 7, units = 'in')
    par(mfrow=c(1,1),oma=c(1,1,1,1), mgp=c(2, 1,0))
    plot(timeseries[,1],  timeseries[,2], type="b", col="black", pch=16, xlab="year", ylim=c(0,max_index*1.2), ylab=dep_text, main=main.lab) # ylim=c(0,max_index*1.2)
    lines(timeseries[,1], (timeseries[,2]-1.96*timeseries[,3]), type="l",lty=2, col="red" )
    lines(timeseries[,1], (timeseries[,2]+1.96*timeseries[,3]), type="l",lty=2, col="red" )
    legend("topright", c("time series", "CI"), lty=c(1,2), pch=c(16, NA), col=c("black","red"))
    dev.off()
} else {
    par(mfrow=c(1,1),oma=c(1,1,1,1), mgp=c(2, 1,0))
    plot(timeseries[,1],  timeseries[,2], type="b", col="black", pch=16, xlab="year", ylim=c(0,max_index*1.2), ylab=dep_text, main=main.lab) # ylim=c(0,max_index*1.2)
    lines(timeseries[,1], (timeseries[,2]-1.96*timeseries[,3]), type="l",lty=2, col="red" )
    lines(timeseries[,1], (timeseries[,2]+1.96*timeseries[,3]), type="l",lty=2, col="red" )
    legend("topright", c("time series", "CI"), lty=c(1,2), pch=c(16, NA), col=c("black","red"))
}
    message("\n Estimation of biomass indices completed ")
    return(list(timeseries_abundance,timeseries_biomass, tab.strata.abu, tab.strata.bio))
}