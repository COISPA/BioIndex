#' Length Frequency Distribution (LFD) Estimation and Plotting
#'
#' This function estimates and plots the length frequency distribution (LFD) by year and by stratum,
#' based on merged biological and haul data from MEDITS surveys. It computes raised numbers,
#' applies stratification weights, and generates publication-ready plots and CSV outputs.
#' Stratification weighting is based on swept area and surface area of the strata.
#'
#' @description Input parameters include biological and haul data, stratification schemes,
#' reference filters (e.g., sex, GSA, country, depth), and output control options (working directory, saving, verbosity).
#'
#' @param mTATC Data frame resulting from the merge of TA and TC datasets using `merge_TATBTC()`. It must include raised numbers.
#' @param sex Character. Target sex for the analysis. Allowed values: `"F"`, `"M"`, `"I"`, `"N"`, or `"all"` (default) for combined sexes.
#' @param GSA Numeric. The GSA (Geographical Sub-Area) code of reference for the analysis.
#' @param country Character vector. Reference country or countries for the analysis. Use `"all"` (default) to include all countries available in the data.
#' @param depth_range Numeric vector of length 2, specifying the minimum and maximum depths (in meters) to filter the hauls used in the analysis.
#' @param strata_scheme Data frame containing the stratification scheme. Must include `CODE`, `MIN_DEPTH`, `MAX_DEPTH`, `GSA`, `COUNTRY`.
#' @param stratification Data frame with surface areas per stratum. Must include columns: `GSA`, `CODE`, `COUNTRY`, `SURF`.
#' @param wd Character. Working directory used to save output plots and tables. Required if `save = TRUE`.
#' @param save Logical. If `TRUE` (default), the function saves output tables and plots to `wd/output/`. If `FALSE`, no files are saved.
#' @param verbose Logical. If `TRUE` (default), informative messages are printed to the console to trace the function steps.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{`LFD`}{A data frame of total LFD by year across all strata.}
#'   \item{`LFD by stratum`}{A data frame of LFD by year and by stratum.}
#' }
#'
#' @details
#' The function:
#' \itemize{
#'   \item Filters hauls by depth and country.
#'   \item Computes raised numbers per haul.
#'   \item Applies stratified weights based on the area per stratum.
#'   \item Outputs and/or plots the length frequency distribution by year and stratum.
#'   \item Handles sex-specific or combined-sex analyses.
#' }
#'
#' The plot outputs include:
#' \itemize{
#'   \item Combined LFD by year across all strata.
#'   \item LFD by stratum (faceted), stratified by year.
#' }
#'
#' Output files are saved as CSV and high-resolution JPEGs if `save = TRUE` and `wd` is defined.
#'
#' @importFrom dplyr summarise summarize group_by arrange
#' @importFrom reshape2 melt dcast
#' @importFrom ggplot2 ggplot geom_line xlab ylab ggtitle xlim ylim theme_bw aes facet_wrap
#' @export
LFD <- function(mTATC, sex="all", GSA, country="all", depth_range, strata_scheme, stratification, wd=NA, save=TRUE, verbose=TRUE){
    if (is.na(wd)) { wd <- tempdir() }

    id <- NULL
    if (FALSE) {

        wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\R_BioIndex_3.3_(in update)"
        # ta <- read.table(paste(wd,"input/TA/TA GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
        # tb <- read.table(paste(wd,"input/TB/TB GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
        # tc <- read.table(paste(wd,"input/TC/TC GSA18 2017-2020.csv", sep="/"), sep=";", header=T)

        ta <- read.table(paste("D:\\OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L\\B-USEFUL\\_MEDITS_Corrected_data_\\TA_GSA18_1999_2021_completo.csv", sep="/"), sep=";", header=T)
        tb <- read.table(paste("D:\\OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L\\B-USEFUL\\_MEDITS_Corrected_data_\\TB_GSA18_1999_2021_completo.csv", sep="/"), sep=";", header=T)
        tc <- read.table(paste("D:\\OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L\\B-USEFUL\\_MEDITS_Corrected_data_\\TC_GSA18_1999_2021_completo.csv", sep="/"), sep=";", header=T)
        # ta <- ta[ta$YEAR==2017, ]
        # tb <- tb[tb$YEAR==2017, ]
        # tc <- tc[tc$YEAR==2017, ]
        depth_range=c(10,800)
        country="all"
        sex <- "all"
        species <- "MERLMER"
        m <- merge_TATBTC(ta, tb, tc, species=species, country=country, wd=wd, verbose=TRUE)
        mTATC <- m[[2]]

        LFD(mTATC, sex=sex, GSA=18, country=country, depth_range=c(10,800), strata_scheme, stratification, wd, save=TRUE)

    }


    if (is.na(wd) & save) {
        save=TRUE
        if (verbose){
            message("Missing working directory. Results are not saved in the local folder.")
        }
    }

    Abundance <- LENGTH_CLASS <- STRATUM <- STRATUM_CODE <- SWEPT_AREA <- YEAR <-
    analysis_stratum1 <- analysis_stratum2 <- analysis_stratum3 <-
    analysis_stratum4 <- analysis_stratum5 <- analysis_stratum6 <- n_raised <-
    num_weighted <- NULL

    mTATC$n_raised <- mTATC$NUMBER_OF_INDIVIDUALS_IN_THE_LENGTH_CLASS_AND_MATURITY_STAGE * mTATC$RAISING_FACTOR

    ddd <-  mTATC
    class_code <- unique(ddd[ddd$GENUS != -1, "LENGTH_CLASSES_CODE"])
    if (class_code == "m") {LC <- 1} else {
        if (class_code == "0") {LC <- 5} else {
            if (class_code == "1") {LC <- 10}
        }
    }


    GENERE <- as.character(unique(mTATC$GENUS)[unique(mTATC$GENUS) != -1])
    SPECIE <- as.character(unique(mTATC$SPECIES)[unique(mTATC$SPECIES) != -1])
    sspp <- paste(GENERE,SPECIE, sep="")
    GSA <- unique(mTATC$GSA)
    species <- sspp



    countries <- unique(mTATC[!is.na(mTATC$COUNTRY),"COUNTRY"])
    l_country <- length(unique(mTATC$COUNTRY))

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

    strata_scheme <- strata_scheme[strata_scheme$GSA == GSA & strata_scheme$COUNTRY %in% as.character(unique(country_analysis))[1], ]

    stratification$SURF <- as.numeric(stratification$SURF)

    depth <- data.frame(matrix(NA,ncol=4, nrow=length(strata_scheme$CODE)))
    colnames(depth) <- c("strata", "min", "max", "bul")
    depth$strata <- strata_scheme$CODE
    depth$min <-strata_scheme$MIN_DEPTH
    depth$max <-strata_scheme$MAX_DEPTH
    if (depth_range[2] != 800) {depth_range[2] <- depth_range[2]}
    mTATC <-  mTATC[mTATC$MEAN_DEPTH>depth_range[1] & mTATC$MEAN_DEPTH<=depth_range[2],]
    ddd <-  ddd[ddd$MEAN_DEPTH>depth_range[1] & ddd$MEAN_DEPTH<=depth_range[2],]
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
    strata_scheme$WEIGHT <- 0
    if ("1" %in% strata_scheme$CODE) {strata_scheme[strata_scheme$CODE =="1", "WEIGHT"] <- peso_s1}
    if ("2" %in% strata_scheme$CODE) {strata_scheme[strata_scheme$CODE =="2", "WEIGHT"] <- peso_s2}
    if ("3" %in% strata_scheme$CODE) {strata_scheme[strata_scheme$CODE =="3", "WEIGHT"] <- peso_s3}
    if ("4" %in% strata_scheme$CODE) {strata_scheme[strata_scheme$CODE =="4", "WEIGHT"] <- peso_s4}
    if ("5" %in% strata_scheme$CODE) {strata_scheme[strata_scheme$CODE =="5", "WEIGHT"] <- peso_s5}
    if ("6" %in% strata_scheme$CODE) {strata_scheme[strata_scheme$CODE =="6", "WEIGHT"] <- peso_s6}


    ddd <- mTATC[mTATC$LENGTH_CLASS != -1,]
    if (!sex %in% c("I","N","M","F","all")) {
        stop("Unexpected sex selection for the analysis.\n")
    } else if (sex=="I") {
        ddd <- ddd[ddd$SEX == "I", ]
        # mTATC <- mTATC[mTATC$SEX == "I", ]
    } else if (sex=="N") {
        ddd <- ddd[ddd$SEX == "N", ]
        # mTATC <- mTATC[mTATC$SEX == "N", ]
    } else if (sex=="M") {
        ddd <- ddd[ddd$SEX == "M", ]
        # mTATC <- mTATC[mTATC$SEX == "M", ]
    } else if (sex=="F") {
        ddd <- ddd[ddd$SEX == "F", ]
        # mTATC <- mTATC[mTATC$SEX == "F", ]
    }

    if (nrow(ddd) != 0) {
        ddd <- ddd[ddd$GSA %in% GSA & ddd$COUNTRY %in% as.character(country_analysis), ]
        mTATC <- mTATC[mTATC$GSA %in% GSA & mTATC$COUNTRY %in% as.character(country_analysis), ]
    }

    if (nrow(ddd) != 0 & nrow(mTATC) != 0) {

        df1 <- suppressMessages( data.frame(ddd %>% group_by(YEAR, STRATUM_CODE,LENGTH_CLASS) %>% summarise(numb_per_stratum = sum(n_raised))))
        df1$id <- paste(df1$YEAR,df1$STRATUM_CODE,sep="_")

        df2t <-  suppressMessages(data.frame(mTATC %>% group_by(id, YEAR, STRATUM_CODE) %>% summarise(SWEPT_AREA=mean(SWEPT_AREA))))
        df2t$id <- paste(df2t$YEAR,df2t$STRATUM_CODE,sep="_")
        df2 <- suppressMessages(data.frame(df2t %>% group_by(id, YEAR, STRATUM_CODE) %>% summarise(SWEPT_AREA=sum(SWEPT_AREA))))

        mdf <- merge(df1,df2[, c("id", "SWEPT_AREA")], by.x="id",by.y="id",all=TRUE)
        mdf$n_km <- mdf$numb_per_stratum / mdf$SWEPT_AREA
        mw <- merge(mdf,strata_scheme[, c("CODE","WEIGHT")],by.x="STRATUM_CODE", by.y="CODE",all.x=TRUE)
        mw$num_weighted <- mw$n_km * mw$WEIGHT
        mw <- mw[!is.na(mw$STRATUM_CODE),]


        df <- dcast(mw, LENGTH_CLASS + YEAR ~ STRATUM_CODE, value.var = "n_km")
        strata <- colnames(df[3:ncol(df)])
        strata <- paste("STRATUM",strata,sep=" ")
        colnames(df) <- c("LC","YEAR", strata)
        LC_range <- range(mw$LENGTH_CLASS)
        LCs <- seq(LC_range[1],LC_range[2],LC)
        LCs <- data.frame(LC = LCs)
        years <- sort(unique(df$YEAR))
        msc <- list()
        y=1
        for (y in 1:length(years)) {
            temp <- merge(LCs,df[df$YEAR == years[y], ], by.x="LC",by.y="LC",all.x=TRUE, all.y=TRUE)
            temp$YEAR <- years[[y]]
            temp[is.na(temp)] <- 0
            msc[[y]] <- temp
        }

        mst <- do.call(rbind,msc)
        mst$id <- paste(mst$YEAR,mst$LC,sep="_")
        df_all <- suppressMessages(data.frame(mw %>% group_by(YEAR,LENGTH_CLASS) %>% summarise(`ALL STRATA` = sum(num_weighted))))
        df_all$id <- paste(df_all$YEAR,df_all$LENGTH_CLASS,sep="_")
        m <- merge(mst,df_all[, c("ALL.STRATA","id")], by.x="id",by.y="id", all.x=TRUE, all.y=TRUE)
        m[is.na(m)] <- 0

        m2 <- m[, -c(which(colnames(m)=="id"))]
        m2 <- m2 %>% arrange(YEAR,LC)
        m2all <- m2[, c("LC","YEAR","ALL.STRATA")]
        m2all <- dcast(m2all, LC ~ YEAR, value.var = "ALL.STRATA")

        m2s <-  m2[, -(which(colnames(m2)=="ALL.STRATA"))]
        m2s <- melt(m2s, id.vars=c("LC","YEAR"),variable.name = "STRATUM",value.name = "Abundance")
        m2s <- dcast(m2s, LC + STRATUM ~ YEAR, value.var = "Abundance")
        m2s <- m2s %>% arrange(STRATUM, LC)
        ss <- strsplit(as.character(m2s$STRATUM)," ")
        ss <- do.call(rbind,ss)[,2]
        m2s$STRATUM <- ss

        if (verbose) message("LFD successfully estimated.")

     if (save) {
        if (sex=="all"){
        write.table(m2s,paste(wd, "/output/",sspp,"_GSA",GSA,"_LFD_(Combined_by_stratum).csv", sep=""),sep=";",row.names=F)
        write.table(m2all,paste(wd, "/output/",sspp,"_GSA",GSA,"_LFD_(Combined).csv", sep=""),sep=";",row.names=F)
        if (verbose) message("LFD tables saved in output folder.")
        } else {
        write.table(m2s,paste(wd, "/output/",sspp,"_GSA",GSA,"_LFD_(",sex,"_by_stratum).csv", sep=""),sep=";",row.names=F)
        write.table(m2all,paste(wd, "/output/",sspp,"_GSA",GSA,"_LFD_(",sex,").csv", sep=""),sep=";",row.names=F)
        if (verbose) message("LFD tables saved in output folder.")
            }
     }

        m <- melt(m, id.vars=c("id","LC","YEAR"),variable.name = "STRATUM",value.name = "Abundance")
        m$YEAR <- as.factor(m$YEAR)

# PLOT by year
        if (sex=="all"){
            main <- paste(sspp,"_GSA",GSA,"_(LFD_by_year_Combained)_",depth_range[1],"-",depth_range[2], " m", sep="")
            main.lab <- paste(sspp," GSA",GSA," (LFD by year, Combained) ",depth_range[1],"-",depth_range[2], " m", sep="")
            dep_text <-expression(paste("LFD", sep=" "))
            min_lc <- min(m[m$STRATUM == "ALL.STRATA", "LC"])
            max_lc <- max(m[m$STRATUM == "ALL.STRATA", "LC"])
            max_freq <- max(m[m$STRATUM == "ALL.STRATA", "Abundance"])

            p <- ggplot(data=m[m$STRATUM == "ALL.STRATA", ], aes(x=LC, y=Abundance, group=YEAR, colour=YEAR)) +
                geom_line() +
                xlab("LC") + ylab(expression(paste("Abundance ", (n/km^2), sep=" "))) +
                ggtitle(main.lab) +
                xlim(min_lc, max_lc) +
                ylim(0,max_freq) +
                theme_bw()
            print(p)
            if(save){
                if(save) ggsave(paste(wd,"/output/",main,".jpg",sep=""), dpi=300 , width=6, height=5, plot = p) #
                }
        } else {
            main <- paste(sspp,"_GSA",GSA,"_(LFD_by_year_",toupper(sex),")_",depth_range[1],"-",depth_range[2], " m", sep="")
            main.lab <- paste(sspp," GSA",GSA," (LFD by year, ",toupper(sex),") ",depth_range[1],"-",depth_range[2], " m", sep="")
            dep_text <-expression(paste("LFD", sep=" "))
            min_lc <- min(m[m$STRATUM == "ALL.STRATA", "LC"])
            max_lc <- max(m[m$STRATUM == "ALL.STRATA", "LC"])
            max_freq <- max(m[m$STRATUM == "ALL.STRATA", "Abundance"])

            p <- ggplot(data=m[m$STRATUM == "ALL.STRATA", ], aes(x=LC, y=Abundance, group=YEAR, colour=YEAR)) +
                geom_line() +
                xlab("LC") + ylab(expression(paste("Abundance ", (n/km^2), sep=" "))) +
                ggtitle(main.lab) +
                xlim(min_lc, max_lc) +
                ylim(0,max_freq) +
                theme_bw()
            print(p)
            if (save){
                if(save) ggsave(paste(wd,"/output/",main,".jpg",sep=""), dpi=300 , width=6, height=5, plot = p)
                }
        }



        # plot by strata

        if (sex=="all"){
            main <- paste(sspp,"_GSA",GSA,"_(LFD_by_stratum_Combained)_",depth_range[1],"-",depth_range[2], " m", sep="")
            main.lab <- paste(sspp," GSA",GSA," (LFD by stratum, Combained) ",depth_range[1],"-",depth_range[2], " m", sep="")
            dep_text <-expression(paste("LFD", sep=" "))
            min_lc <- min(m[m$STRATUM != "ALL.STRATA", "LC"])
            max_lc <- max(m[m$STRATUM != "ALL.STRATA", "LC"])
            max_freq <- max(m[m$STRATUM != "ALL.STRATA", "Abundance"])

            p2 <- ggplot(data=m[m$STRATUM != "ALL.STRATA", ], aes(x=LC, y=Abundance, group=YEAR, colour=YEAR)) +
                geom_line(linewidth=0.3) +
                xlab("LC") + ylab(expression(paste("Abundance ", (n/km^2), sep=" "))) +
                ggtitle(main.lab) +
                xlim(min_lc, max_lc) +
                ylim(0,max_freq) +
                theme_bw()+
                facet_wrap(~ STRATUM)  # , ncol=cols.graph
            print(p2)
            if (save) {
                if(save) ggsave(paste(wd,"/output/",main,".jpg",sep=""), dpi=300 , width=6, height=5, plot = p2) #
            }
        } else {
            main <- paste(sspp,"_GSA",GSA,"_(LFD_by_stratum_",toupper(sex),")_",depth_range[1],"-",depth_range[2], " m", sep="")
            main.lab <- paste(sspp," GSA",GSA," (LFD by stratum, ",toupper(sex),") ",depth_range[1],"-",depth_range[2], " m", sep="")
            dep_text <-expression(paste("LFD", sep=" "))
            min_lc <- min(m[m$STRATUM != "ALL.STRATA", "LC"])
            max_lc <- max(m[m$STRATUM != "ALL.STRATA", "LC"])
            max_freq <- max(m[m$STRATUM != "ALL.STRATA", "Abundance"])

            p2 <- ggplot(data=m[m$STRATUM != "ALL.STRATA", ], aes(x=LC, y=Abundance, group=YEAR, colour=YEAR)) +
                geom_line(linewidth=0.3) +
                xlab("LC") + ylab(expression(paste("Abundance ", (n/km^2), sep=" "))) +
                ggtitle(main.lab) +
                xlim(min_lc, max_lc) +
                ylim(0,max_freq) +
                theme_bw()+
                facet_wrap(~ STRATUM)  # , ncol=cols.graph
            print(p2)
            if(save){
                if(save) ggsave(paste(wd,"/output/",main,".jpg",sep=""), dpi=300 , width=6, height=5, plot = p2)
            }
        }

return(list("LFD"=m2all,"LFD by stratum"=m2s))
    } else {
        if (verbose) message("Not enough data for LFD estimation.")
    }

}
