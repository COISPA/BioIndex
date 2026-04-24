#' Spearman test
#'
#' @param abundance data frame of abundance indices
#' @param biomass data frame of biomass indices
#' @param years reference years for the analysis
#' @param sspp reference species for the analysis
#' @param wd path of working directory
#' @param save boolean. If TRUE the plot is saved in the user defined working directory (wd)
#' @importFrom utils write.table
#' @export
spearman <- function(abundance=NA, biomass=NA, years, sspp=NA, wd=NA, save=TRUE){

    if (FALSE) {

        # library(BioIndex)
        years <- c(2017, 2020)

        GSA=18
        save=TRUE
        country="all"
        depth_range <- c(10,800)
        species <- "MERLMER"

        wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\R_BioIndex_3.3_(in update)"
        verbose=TRUE
        wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\R_BioIndex_3.3_(in update)"
        ta <- read.table(paste(wd,"input/TA/TA GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
        tb <- read.table(paste(wd,"input/TB/TB GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
        tc <- read.table(paste(wd,"input/TC/TC GSA18 2017-2020.csv", sep="/"), sep=";", header=T)


        m <- merge_TATBTC(ta, tb, tc, species=species, country=country, wd=wd, verbose=TRUE)
        mTATB <- m[[1]]
        index <- indices_ts(mTATB, GSA=GSA, country=country, depth_range=depth_range, strata_scheme=BioIndex::strata_scheme, stratification=BioIndex::stratification,wd, save=FALSE)

        abundance= index[[1]]
        biomass= index[[2]]

        spearman(abundance=abundance, biomass=biomass, years, sspp=species, wd=wd, save=TRUE)
    }



    ########   START   ########

    if (is.na(wd) & save) {
        save =FALSE
        if (verbose){
            message("Missing working directory. Results are not saved in the local folder.")
        }
    }


    yrange <- data.frame(range = as.numeric(range(years)))
    rownames(yrange) <- c("start", "end")


#----------------------------
# Abundance indices
#----------------------------
    if(all(!is.na(abundance)) & class(abundance)%in% c("data.frame")){
        timeseries <- abundance
        years_sel <- (yrange[1,1]:yrange[2,1])
        timeseries <- timeseries[timeseries$year %in% years_sel , ]
        x=log(timeseries[,2]+1)
        spear(x)
        spearman_abundance <- spear(x)
        spearman_abundance <- data.frame(index="abundance", spearman_abundance)
        if(is.infinite(spearman_abundance$t)) {spearman_abundance$t <- "-Inf"}
    } else {
        spearman_abundance <- data.frame(matrix(ncol=4,nrow=1))
        colnames(spearman_abundance) <- c("index","r","t","p")
        spearman_abundance[1,] <- c("abundance",NA,NA,NA)

    }

#----------------------------
# Biomass indices
#----------------------------
    if(all(!is.na(biomass)) & class(biomass)%in% c("data.frame")){
        timeseries <- biomass
        years_sel <- (yrange[1,1]:yrange[2,1])
        timeseries <- timeseries[timeseries$year %in% years_sel , ]
        x=log(timeseries[,2]+1)
        spear(x)
        spearman_biomass <- spear(x)
        spearman_biomass <- data.frame(index="biomass", spearman_biomass)
        if(is.infinite(spearman_biomass$t)) {spearman_biomass$t <- "-Inf"}
    } else {
        spearman_biomass <- data.frame(matrix(ncol=4,nrow=1))
        colnames(spearman_biomass) <- c("index","r","t","p")
        spearman_biomass[1,] <- c("biomass",NA,NA,NA)
    }
#----------------------------
# Results in table
#----------------------------
    tab_spearman <- rbind(spearman_abundance , spearman_biomass)

    if(save & !is.na(wd)){
        write.table(tab_spearman, paste(wd, "/output/",sspp," - Spearman summary_RSS.csv", sep=""), sep=";", row.names=FALSE)
    }

    return(tab_spearman)

}
