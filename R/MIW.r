#' Estimation of Mean Individual Weight (MIW) time series
#'
#' @param mTATB data frame of the merged TA and TB
#' @param GSA reference GSA for the analysis
#' @param country reference countries in the GSA for the analysis
#' @param depth_range range of depth strata to perform the analysis (min, max)
#' @param strata_scheme data frame of the stratification scheme
#' @param stratification data frame of strata surface area
#' @param wd working directory
#' @param save boolean. If TRUE the plot is saved in the user defined working directory (wd)
#' @param verbose boolean. If TRUE messages are reported in the console
#' @importFrom gridExtra grid.arrange
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_density ggtitle geom_boxplot xlab ylab ylim xlim ggsave theme_bw geom_line geom_point
#' @importFrom stringr str_split
#' @import grDevices
#' @import graphics
#' @import magrittr
#' @export
#'
#'
MIW <- function(mTATB, GSA, country="all", depth_range, strata_scheme, stratification, wd=NA, save=TRUE,verbose=TRUE){

    if (FALSE) {
        GSA <- 18
        country <- "all"
        depthr <- c(10,800)
        # strata=BioIndex::strata_scheme
        # stratification_tab = BioIndex::stratification
        strata_scheme=BioIndex::strata_scheme
        stratification = BioIndex::stratification
        # wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\Test_BioIndex_package"
        save=FALSE
        verbose=TRUE

        # MIW(mTATB, GSA, country="all", depth_range=depthr, strata_scheme=strata, stratification=stratification_tab, wd, save,verbose)
    }

#-----------------------
# MIW INDICES
#-----------------------
year <- value <- NULL
#---- COUNTRIES
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
merge_TATB <- mTATB
merge_TATB <- merge_TATB[merge_TATB$COUNTRY %in% country_analysis, ]

strata_scheme <- strata_scheme[strata_scheme$GSA == GSA & strata_scheme$COUNTRY %in% as.character(unique(country_analysis))[1], ]

#---
dependent <- "MIW"
dep_text <-"MIW (kg)"
GENERE <- as.character(unique(merge_TATB$GENUS)[unique(merge_TATB$GENUS) != -1])
SPECIE <- as.character(unique(merge_TATB$SPECIES)[unique(merge_TATB$SPECIES) != -1])
sspp <- paste(GENERE,SPECIE, sep="")
GSA <- unique(merge_TATB$GSA)
species <- sspp
#---
data <- merge_TATB
depth <- data.frame(matrix(NA,ncol=4, nrow=length(strata_scheme$CODE)))
colnames(depth) <- c("strata", "min", "max", "bul")
depth$strata <- strata_scheme$CODE # c(1,2,3,4,5)
depth$min <-strata_scheme$MIN_DEPTH  # c(10,50,100,200,500)
depth$max <-strata_scheme$MAX_DEPTH  # c(50,100,200,500,800)
# print("Select the depth range for the analysis")
# depth_range <- dlgInput("Depth range for the analysis: ", default="5,35",Sys.info()[""])$res
# depth_range <- data.frame(range = strsplit(depth_range, ",")); depth_range <- as.numeric(as.character(depth_range[,1]))
if (depth_range[2] != 800) {depth_range[2] <- depth_range[2]}
data <-  data[data$MEAN_DEPTH>depth_range[1] & data$MEAN_DEPTH<=depth_range[2],]
strata_range <- seq(which(depth[,"min"]==depth_range[1]),which(depth[,"max"]==depth_range[2]),1)
depth[depth$strata %in% strata_range, "bul"] <- T
depth[!(depth$strata %in% strata_range), "bul"] <- F
#---
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
        #index computation
        s1 <- s1[!is.na(s1$TOTAL_NUMBER_IN_THE_HAUL ) & s1$TOTAL_NUMBER_IN_THE_HAUL>0, ]
        s1$mean_MIW <- (s1[ , "TOTAL_WEIGHT_IN_THE_HAUL" ]/1000) / s1[ , "TOTAL_NUMBER_IN_THE_HAUL" ]
        s_MIW<- sum(s1[ ,"TOTAL_WEIGHT_IN_THE_HAUL"]/1000)/sum(s1[ ,"TOTAL_NUMBER_IN_THE_HAUL"])
        s_a <- sum(s1$SWEPT_AREA)
        res_table_cala2[i,2] <- s_MIW

        #se computation
        sq <- (s1$mean_MIW - s_MIW)^2
        sqA <- sq*s1$SWEPT_AREA
        sum_sqA <- sum(sqA)/(length(s1[,1])-1)
        f1 <- s_a/area_s1
        se_table2[i,2] <- (((peso_s1)^2 * sum_sqA)/ s_a)* (1-f1)
        sd_strata[i,2] <- sqrt( ((sum_sqA)/ s_a)* (1-f1))
    } else {area_s1 = 0}

    if (analysis_stratum2 == T){
        #index computation
        s2 <- s2[!is.na(s2$TOTAL_NUMBER_IN_THE_HAUL ) & s2$TOTAL_NUMBER_IN_THE_HAUL>0, ]
        s2$mean_MIW <- (s2[ , "TOTAL_WEIGHT_IN_THE_HAUL" ]/1000) / s2[ , "TOTAL_NUMBER_IN_THE_HAUL" ]
        s_MIW<- sum(s2[ ,"TOTAL_WEIGHT_IN_THE_HAUL"]/1000)/sum(s2[ ,"TOTAL_NUMBER_IN_THE_HAUL"])
        s_a <- sum(s2$SWEPT_AREA)
        res_table_cala2[i,3] <- s_MIW

        #se computation
        sq <- (s2$mean_MIW - s_MIW)^2
        sqA <- sq*s2$SWEPT_AREA
        sum_sqA <- sum(sqA)/(length(s2[,1])-1)
        f2 <- s_a/area_s2
        se_table2[i,3] <- (((peso_s2)^2 * sum_sqA)/ s_a)* (1-f2)
        sd_strata[i,3] <- sqrt( ((sum_sqA)/ s_a)* (1-f2))
    } else {area_s2 = 0}

    if (analysis_stratum3 == T){
        #index computation
        s3 <- s3[!is.na(s3$TOTAL_NUMBER_IN_THE_HAUL ) & s3$TOTAL_NUMBER_IN_THE_HAUL>0, ]
        s3$mean_MIW <- (s3[ , "TOTAL_WEIGHT_IN_THE_HAUL" ]/1000) / s3[ , "TOTAL_NUMBER_IN_THE_HAUL" ]
        s_MIW<- sum(s3[ ,"TOTAL_WEIGHT_IN_THE_HAUL"]/1000)/sum(s3[ ,"TOTAL_NUMBER_IN_THE_HAUL"])
        s_a <- sum(s3$SWEPT_AREA)
        res_table_cala2[i,4] <- s_MIW

        #se computation
        sq <- (s3$mean_MIW - s_MIW)^2
        sqA <- sq*s3$SWEPT_AREA
        sum_sqA <- sum(sqA)/(length(s3[,1])-1)
        f3 <- s_a/area_s3
        se_table2[i,4] <- (((peso_s3)^2 * sum_sqA)/ s_a)* (1-f3)
        sd_strata[i,4] <- sqrt( ((sum_sqA)/ s_a)* (1-f3))
    } else {area_s3 = 0}

    if (analysis_stratum4 == T){
        #index computation
        s4 <- s4[!is.na(s4$TOTAL_NUMBER_IN_THE_HAUL ) & s4$TOTAL_NUMBER_IN_THE_HAUL>0, ]
        s4$mean_MIW <- (s4[ , "TOTAL_WEIGHT_IN_THE_HAUL" ]/1000) / s4[ , "TOTAL_NUMBER_IN_THE_HAUL" ]
        s_MIW<- sum(s4[ ,"TOTAL_WEIGHT_IN_THE_HAUL"]/1000)/sum(s4[ ,"TOTAL_NUMBER_IN_THE_HAUL"])
        s_a <- sum(s4$SWEPT_AREA)
        res_table_cala2[i,5] <- s_MIW

        #se computation
        sq <- (s4$mean_MIW - s_MIW)^2
        sqA <- sq*s4$SWEPT_AREA
        sum_sqA <- sum(sqA)/(length(s4[,1])-1)
        f4 <- s_a/area_s4
        se_table2[i,5] <- (((peso_s4)^2 * sum_sqA)/ s_a)* (1-f4)
        sd_strata[i,5] <- sqrt( ((sum_sqA)/ s_a)* (1-f4))
    } else {area_s4 = 0}

    if (analysis_stratum5 == T){
        #index computation
        s5 <- s5[!is.na(s5$TOTAL_NUMBER_IN_THE_HAUL ) & s5$TOTAL_NUMBER_IN_THE_HAUL>0, ]
        s5$mean_MIW <- (s5[ , "TOTAL_WEIGHT_IN_THE_HAUL" ]/1000) / s5[ , "TOTAL_NUMBER_IN_THE_HAUL" ]
        s_MIW<- sum(s5[ ,"TOTAL_WEIGHT_IN_THE_HAUL"]/1000)/sum(s5[ ,"TOTAL_NUMBER_IN_THE_HAUL"])
        s_a <- sum(s5$SWEPT_AREA)
        res_table_cala2[i,6] <- s_MIW

        #se computation
        sq <- (s5$mean_MIW - s_MIW)^2
        sqA <- sq*s5$SWEPT_AREA
        sum_sqA <- sum(sqA)/(length(s5[,1])-1)
        f5 <- s_a/area_s5
        se_table2[i,6] <- (((peso_s5)^2 * sum_sqA)/ s_a)* (1-f5)
        sd_strata[i,6] <- sqrt( ((sum_sqA)/ s_a)* (1-f5))
    } else {area_s5 = 0}

    if (analysis_stratum6 == T){
        #index computation
        s6 <- s6[!is.na(s6$TOTAL_NUMBER_IN_THE_HAUL ) & s6$TOTAL_NUMBER_IN_THE_HAUL>0, ]
        s6$mean_MIW <- (s6[ , "TOTAL_WEIGHT_IN_THE_HAUL" ]/1000) / s6[ , "TOTAL_NUMBER_IN_THE_HAUL" ]
        s_MIW<- sum(s6[ ,"TOTAL_WEIGHT_IN_THE_HAUL"]/1000)/sum(s6[ ,"TOTAL_NUMBER_IN_THE_HAUL"])
        s_a <- sum(s6$SWEPT_AREA)
        res_table_cala2[i,7] <- s_MIW

        #se computation
        sq <- (s6$mean_MIW - s_MIW)^2
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

timeseries <- data.frame(year = res_table_cala2[,1], MIW = res_table_cala2[,8], sd= se_table2[, "sd"])
timeseries$CV <- timeseries[,3] / timeseries[,2]
timeseries$invCV <- 1/timeseries$CV
colnames(res_table_cala2) <- c("year", "stratum 1","stratum 2", "stratum 3", "stratum 4", "stratum 5", "stratum 6", "Indices")



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
    tab.strata.MIW <- data.frame(year=res_table_cala2[,"year"], res_table_cala2[, cols.strata], TAB_SD_CV[,-1])
    colnames(tab.strata.MIW) <- c("year", paste("MIW", depth[!is.na(depth$strata), 'strata'],sep="_"), paste("SD", depth[!is.na(depth$strata), 'strata'],sep="_"), paste("CV", depth[!is.na(depth$strata), 'strata'],sep="_"))
    title.MIW.haul <- paste(sspp,"_GSA",GSA,"_(",dependent,"_by_stratum)-RSS_",depth_range[1],"-",depth_range[2], " m", sep="")


    is.nan.data.frame <- function(x){do.call(cbind, lapply(x, is.nan))}
    tab.strata.MIW[is.nan.data.frame(tab.strata.MIW)] <- NA

    if(save){
    write.table(tab.strata.MIW, paste(wd,"/output/",title.MIW.haul,".csv", sep=""), sep=";", row.names = F)
    }

    sel.strata1 <- c(1,depth[(depth$bul), "strata"] + 1)
    sel.strata2 <- c(1, depth[(depth$bul), "strata"] + 1 + length(cols.strata))
    tab.ind <- tab.strata.MIW[, sel.strata1]
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

    years <- sort(unique(df.plot1$year))
    if(length(years)>1){
         pm1 <- ggplot(data=df.plot1, aes(x=year, y=value, group=strata, colour=strata)) +
            geom_line(na.rm=TRUE) +
            geom_point(na.rm=TRUE) +
            xlab("year") +
            ylab(dep_text) +
            ggtitle(main.lab) +
            theme_bw()
    } else {
         pm1 <- ggplot(data=df.plot1, aes(x=year, y=value, group=strata, colour=strata)) +
            geom_point(na.rm=TRUE) +
            xlab("year") +
            ylab(dep_text) +
            ggtitle(main.lab) +
            theme_bw()
    }

    print(pm1)

    if (save){
        ggsave(paste(wd,"/output/",main,"_Timeseries.jpg",sep=""), dpi=300 , width=6, height=5)
    }


# plot timeseries of mean indices
main <- paste(sspp,"_GSA",GSA,"_(",dependent,")-Random_Stratified_Sampling_",depth_range[1],"-",depth_range[2], " m", sep="")
main.lab <- paste(sspp," GSA",GSA," (",dependent,")-RSS ",depth_range[1],"-",depth_range[2], " m", sep="")
max_index <- max(timeseries[,"MIW"]) + max(timeseries[!is.na(timeseries$sd),"sd"])*1.2
if (save){
write.table(timeseries, paste(wd,"/output/",main,"_Timeseries.csv", sep=""), sep=";", row.names = F)
}

if (save){
    jpeg(paste(wd,"/output/",main,"_Timeseries.jpg",sep=""), res = 300, width = 8, height = 7, units = 'in')
    par(mfrow=c(1,1),oma=c(1,1,1,1), mgp=c(2, 1,0))
    plot(timeseries[,"year"],  timeseries[,"MIW"], type="b", col="black", pch=16, xlab="year", ylim=c(0,max_index*1.2), ylab=dep_text, main=main.lab) # ylim=c(0,max_index*1.2)
    lines(timeseries[,"year"], (timeseries[,"MIW"]-1.96*timeseries[,"sd"]), type="l",lty=2, col="red" )
    lines(timeseries[,"year"], (timeseries[,"MIW"]+1.96*timeseries[,"sd"]), type="l",lty=2, col="red" )
    legend("topright", c("time series", "CI"), lty=c(1,2), pch=c(16, NA), col=c("black","red"))
    dev.off()
} else {
    par(mfrow=c(1,1),oma=c(1,1,1,1), mgp=c(2, 1,0))
    plot(timeseries[,"year"],  timeseries[,"MIW"], type="b", col="black", pch=16, xlab="year", ylim=c(0,max_index*1.2), ylab=dep_text, main=main.lab) # ylim=c(0,max_index*1.2)
    lines(timeseries[,"year"], (timeseries[,"MIW"]-1.96*timeseries[,"sd"]), type="l",lty=2, col="red" )
    lines(timeseries[,"year"], (timeseries[,"MIW"]+1.96*timeseries[,"sd"]), type="l",lty=2, col="red" )
    legend("topright", c("time series", "CI"), lty=c(1,2), pch=c(16, NA), col=c("black","red"))
}
cat("\n Estimation of MIW completed \n")
results <- list(timeseries,tab.strata.MIW)
return(results)
}
