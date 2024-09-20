#' BioIndex
#'
#' @description
#' R code to perform analysis of trawl survey data using MEDITS file format
#'
#' @param ta data frame of the TA table in the MEDITS format
#' @param tb data frame of the TB table in the MEDITS format
#' @param tc data frame of the TC table in the MEDITS format
#' @param sspp reference species for the analysis
#' @param rec_threshold cutoff threshold for recruits
#' @param spaw_threshold cutoff threshold for spawners
#' @param haul_threshold minimum number of individuals to be used in estimation of the spatial indicaticators
#' @param sexes reference sex for the analysis
#' @param depth reference depth range
#' @param GSA reference GSA for the analysis
#' @param country reference country
#' @param map_lim coordinates limits for the maps
#' @param depth_lines depth contours to be plotted in the maps (3 values allowed)
#' @param strata data frame of the reference strata for the study area
#' @param stratification_tab data frame of the stratification scheme
#' @param resolution resolution of the depth line
#' @param buffer buffer around the map
#' @param wd path of the working directory
#' @param zip boolean. If TRUE the results are stored in a zip file into the working directory
#' @param save boolean. If TRUE the results are stored in the working directory
#' @param verbose boolean. If TRUE messages are promted in the console
#' @examples
#' BioIndex(ta=TA, tb=TB, tc=TC, sspp="MERLMER",rec_threshold=200, spaw_threshold=210,sexes="all", depth=c(10,800), GSA=10, country="all", map_lim=c(13.3,15.2,39.9,41.3),depth_lines=c(50,200,800), strata=BioIndex::strata_scheme, stratification_tab = BioIndex::stratification, resolution=1, buffer=0.1, wd=tempdir(), zip=TRUE, save=TRUE, verbose=TRUE)
#'
#' @importFrom methods is
#' @importFrom zip zip
#' @export
BioIndex <- function(ta, tb, tc, sspp,rec_threshold, spaw_threshold,
                     haul_threshold=30, sexes="all", depth, GSA, country="all",
                     map_lim,depth_lines=c(10,200,800),
                     strata=BioIndex::strata_scheme,
                     stratification_tab = BioIndex::stratification,
                     resolution=1, buffer=0.1, wd,zip=TRUE,
                     save=TRUE, verbose=TRUE) {

    if (FALSE) {
        # library(BioIndex)
        wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\Test_BioIndex_package"
        ta <- read.table(file.path(wd,"input","TA GSA18 2017-2020.csv"),sep=";",header=TRUE)
        tb <- read.table(file.path(wd,"input","TB GSA18 2017-2020.csv"),sep=";",header=TRUE)
        tc <- read.table(file.path(wd,"input","TC GSA18 2017-2020.csv"),sep=";",header=TRUE)

        ta <- read.table(file.path(wd,"input","TA.csv"),sep=";",header=TRUE)
        tb <- read.table(file.path(wd,"input","TB.csv"),sep=";",header=TRUE)
        tc <- read.table(file.path(wd,"input","TC.csv"),sep=";",header=TRUE)
        sspp <- "MERLMER"
        rec_threshold=200
        spaw_threshold=210
        haul_threshold=30
        sexes <- "all"
        depth <- c(10,800)
        GSA <- 20
        country <- "all"
        map_lim <- c(15.5,20.0,39.8,42.5)
        depth_lines <- c(200,500,800)
        buffer=0.1
        res=1
        strata=BioIndex::strata_scheme
        stratification_tab = BioIndex::stratification
        save=TRUE
        verbose=TRUE

        BioIndex(ta, tb, tc, sspp,rec_threshold=rec_threshold, spaw_threshold=spaw_threshold,sexes="all", depth=depth, GSA=GSA, country="all", map_lim=map_lim,depth_lines=depth_lines, strata=strata, stratification_tab = stratification_tab, resolution=res, buffer=buffer, wd=wd, zip=TRUE, save=TRUE, verbose=TRUE)


        # ta <- read.table(file.path(wd,"input","TA_Combine.csv"),sep=";",header=TRUE)
        # tb <- read.table(file.path(wd,"input","TB_Combine.csv"),sep=";",header=TRUE)
        # tc <- read.table(file.path(wd,"input","TC_Combine.csv"),sep=";",header=TRUE)
        # sspp <- species <- "PSETMAX"
        # rec_threshold=201
        # spaw_threshold=367
        # haul_threshold=30
        # sexes <- "all"
        # depth <- c(10,125)
        # GSA <- 29
        # country <- "TUR"
        # map_lim <- c(26.9,42.28,40.33,46.89)
        # depth_lines <- c(20,50,100)
        # buffer=0.1
        # res=1
        # stratas=read.table("D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\Test_BioIndex_package\\strata_BS.csv",sep=";",header=TRUE)
        # stratification_tabs = read.table("D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\Test_BioIndex_package\\stratification scheme_BS.csv",sep=";",header=TRUE)
        # strata=stratas
        # stratification = stratification_tabs
        # depth_range=depth
        # save=TRUE
        # verbose=TRUE



        # BioIndex(ta, tb, tc, sspp,rec_threshold=rec_threshold, spaw_threshold=spaw_threshold,sexes="all", depth=depth, GSA=GSA, country="all", map_lim=map_lim,depth_lines=depth_lines, strata=stratas, stratification_tab = stratification_tabs, resolution=res, buffer=buffer, wd=wd, zip=FALSE, save=TRUE, verbose=TRUE)
    }


if (dir.exists(file.path(wd,"output"))){
   dir <- list.files(file.path(wd,"output"), recursive=TRUE, full.names = TRUE, include.dirs = TRUE)
   unlink(dir , recursive=T)
   }



    m <- merge_TATBTC(ta, tb, tc, species=sspp, country=country, wd=wd, strata=strata, verbose=TRUE)

    mTATB <- m[[1]]
    mTATC <- m[[2]]

    ms <- overlayGrid(mTATB=mTATB, mTATC=mTATC, GSA=GSA, country=country, wd=wd, save=save, verbose=verbose)

    mTATBsp <- ms[[1]]
    mTATCsp <- ms[[2]]

    hauls_position(mTATB=mTATB, country=country ,map_lim=map_lim,depth_lines=depth_lines, buffer=buffer, res=resolution, wd=wd,save, verbose)


    #--------------------------------------------------------------
    # Map of abundance and biomass indices by haul
    #--------------------------------------------------------------
    cat("\n########################\n")
    cat("Plot of indices by haul\n")
    cat("########################\n")
    cat("\n")

    source.check <- try(
        bubble_plot_by_haul_indexes(mTATB, map_lim,  depth_lines, buffer=0, res=resolution,wd=wd,save, verbose)
        , silent=T)

    if ( !is(source.check,"try-error") )  {
        cat("\nPlot of indices by haul - completed\n")
    } else {
        message("\nPlot of indices by haul skipped - (Run Error)\n")
    }


    #--------------------------------------------------------------
    # Abundance and biomass indices per GSA in the timeseries
    #--------------------------------------------------------------

    cat("\n########################\n")
    cat("Time series of indices\n")
    cat("########################\n")
    cat("\n")
    index <- indices_ts(mTATB, GSA=GSA, country=country, depth_range=depth, strata_scheme=strata, stratification=stratification_tab,wd, save)
    cat("Time series of indices - completed\n")




    #--------------------------------------------------------------
    # Mean Individual Weights (MIW) per GSA in the timeseries
    #--------------------------------------------------------------

    cat("\n########################\n")
    cat("Time series of indices\n")
    cat("########################\n")
    cat("\n")
    MIW(mTATB, GSA, country, depth_range=depth, strata_scheme=strata, stratification=stratification_tab, wd, save,verbose)
    cat("Time series of MIW - completed\n")



    #--------------------------------------------------------------
    # Sex ratio per GSA in the timeseries
    #--------------------------------------------------------------
    cat("\n########################\n")
    cat("Sex-ratio time series\n")
    cat("########################\n")
    cat("\n")

    SR_analysis <- "ok"

    if (SR_analysis=="ok"){
        if (sum(mTATB$NB_OF_FEMALES + mTATB$NB_OF_MALES, na.rm = TRUE)==0) {
            message("Not enough sex data for sex-ratio estimation")
            cat("\nSex-ratio analysis skipped\n")
        } else {
            sex_ratio(mTATB, GSA, country=country, depth_range=depth, stratas=strata, stratification=stratification_tab, wd, save)
            cat("\nSex-ratio analysis - completed\n")
        }
    } else {
        cat("\nSex-ratio analysis skipped\n")
    }



    #--------------------------------------------------------------
    # Abundance indices of spawners in the time series
    #--------------------------------------------------------------
    cat("\n############################\n")
    cat("Spawners' abundance indices\n")
    cat("############################\n")
    cat("\n")
    #------> check the threshold in the file "~/input/maturity_sizes.csv"
    df_cutoff <- spaw_threshold
    if (is.na(df_cutoff)) {
        message("The SPAWNERS' threshold value not provided. Please, define a value in the 'spaw_threshold' parameter.")
    }
    skip_spawners <- FALSE
    spaw_analysis <-"ok"

    if (spaw_analysis=="ok"){


        source.check <- try(
            index_spawn(mTATB,mTATC, GSA, country, depth_range=depth, cutoff=spaw_threshold, stratification=stratification_tab, wd, save)
            , silent=T)

        if ( !is(source.check,"try-error") )  {
            index_spawn(mTATB,mTATC, GSA, country, depth_range=depth, cutoff=spaw_threshold, stratification=stratification_tab, wd, save)
            if (skip_spawners){
                cat("\nSpawners' indices analysis skipped\n")
            } else {
                cat("\nSpawners' indices analysis - completed\n")
            }

        } else {
            message("\nSpawners' indices analysis skipped - (Run Error)\n")
        }

    }  else {
        cat("\nSpawners' indices analysis skipped\n")
    }


    #--------------------------------------------------------------
    # Abundance indices of recruits in the time series
    #--------------------------------------------------------------
    cat("\n############################\n")
    cat("Recruits' abundance indices\n")
    cat("############################\n")
    cat("\n")
    #------> check the threshold in the file "~/input/maturity_sizes.csv"
    df_cutoff <- rec_threshold
    if (is.na(df_cutoff)) {
        message("The RECRUITS' threshold value not provided. Please, define a value in the 'rec_threshold' parameter.")
    }
    skip_recruits <- FALSE
    rec_analysis <-"ok"

    if (rec_analysis=="ok"){


        source.check <- try(
            index_recr(mTATB,mTATC, GSA, country, depth_range=depth, cutoff=rec_threshold, stratification=stratification_tab, wd, save)
            , silent=T)

        if ( !is(source.check,"try-error") )  {
            index_recr(mTATB,mTATC, GSA, country, depth_range=depth, cutoff=rec_threshold, stratification=stratification_tab, wd, save)
            if (skip_recruits){
                cat("\nRecruits' indices analysis skipped\n")
            } else {
                cat("\nRecruits' indices analysis - completed\n")
            }

        } else {
            message("\nRecruits' indices analysis skipped - (Run Error)\n")
        }

    }  else {
        cat("\nRecruits' indices analysis skipped\n")
    }



    #--------------------------------------------------------------
    # LFD & L0.95
    #--------------------------------------------------------------
    cat("\n############################\n")
    cat("LFD, L0.50 & L0.95\n")
    cat("############################\n")
    cat("\n")

    LFD_analysis <- "ok"

    if (LFD_analysis=="ok"){

        if (sexes == "all") {
            lfd <- LFD(mTATC, sex="all", GSA=18, country=country, depth_range=depth, strata_scheme=strata, stratification=stratification_tab, wd, save)
            Lquant(lfd[[1]], wd, sspp, GSA, save, verbose)
        }

        if (sexes=="M") {
            lfd <- LFD(mTATC, sex="M", GSA=18, country=country, depth_range=depth, strata_scheme=strata, stratification=stratification_tab, wd, save)
            Lquant(lfd[[1]], wd, sspp, GSA, save, verbose)
        }
        if (sexes=="F") {
            lfd <- LFD(mTATC, sex="F", GSA=18, country=country, depth_range=depth, strata_scheme=strata, stratification=stratification_tab, wd, save)
            Lquant(lfd[[1]], wd, sspp, GSA, save, verbose)
        }
        cat("\nLFD & L0.95 analysis - completed\n")
    } else {
        cat("\nLFD & L0.95 analysis skipped\n")
    }


    #--------------------------------------------------------------
    # Spearman test of trends on short timeseries
    #--------------------------------------------------------------
    cat("\n###########################################\n")
    cat("Spearman test of trends on short timeseries\n")
    cat("###########################################\n")
    cat("\n")
    years <- sort(unique(mTATB$YEAR))
    if (length(years) >= 3) {
        Trend_analysis <- "ok"
        if (Trend_analysis=="ok"){
            spearman(abundance=index[[1]], biomass=index[[2]], years, sspp, wd, save)
            cat("\nSpearman test - completed\n")
        } else {
            cat("\nSpearman test skipped\n")
        }

    } else {
        cat("\nSpearman test skipped\n")
    }


    #--------------------------------------------------------------
    # Abundance indices for statistical squares
    # inverse of CV of abundance indices for statistical squares
    # Biomass indices for statistical squares
    # Mean individual weight for statistical squares
    #--------------------------------------------------------------
    depth_range <- paste(depth,collapse = ",")

    index_on_grid(mTATBsp, stratum=depth_range, wd, map_range=map_lim,
                  threshold = haul_threshold,
                  verbose = TRUE, save = TRUE)



    #--------------------------------------------------------------
    # Sex ratio for statistical squares
    #--------------------------------------------------------------
    cat("\n#############################\n")
    cat("Sex-ratio on GFCM grid\n")
    cat("#############################\n")
    cat("\n")

    #------> select the minimum number of individuals per haul to be considered in the analysis
    sexratio_grid_skip <- FALSE
    SR_GRID_analysis <-"ok"
    depth_stratum <- paste(depth,collapse=",")
    if (SR_GRID_analysis=="ok"){
        sex_ratio_on_grid(mTATBsp=mTATBsp, depth=depth_stratum,
                           wd=wd, map_range=map_lim,
                           threshold=haul_threshold,
                           verbose=TRUE,
                           save=TRUE)

        if (sexratio_grid_skip){
            cat("\nSex-ratio on GFCM grid skipped\n")
        } else {
            cat("\nSex-ratio on GFCM grid - completed\n")
        }

    } else {
        cat("\nSex-ratio on GFCM grid skipped\n")
    }




    #--------------------------------------------------------------
    # Bubble plots - indices of recruits (abundance)
    # Bubble plots - indices of spawners (abundance)
    #--------------------------------------------------------------
    cat("\n################################################\n")
    cat("Bubble plots - indices of recruits and spawners\n")
    cat("################################################\n")
    cat("\n")
    #------> check the threshold in the file "~/input/maturity_sizes.csv"

    if (any(is.na(c(rec_threshold, spaw_threshold)))) {
        message("Missing threshold value for the selected species.")
    } else {

    bubbleplot_RS_by_hauls(mTATC=mTATC,
                           map_range=map_lim,
                           thresh_rec=rec_threshold,
                           thresh_spaw=spaw_threshold,
                           depths = depth_lines,
                           buffer=buffer,
                           res=resolution,
                           wd=wd,
                           save = save,
                           verbose = verbose)
        cat("\nBubble plots - indices of recruits and spawners - completed\n")

}








    cat("\n###################\n")
    cat(" Analysis completed\n")
    cat("###################\\n")
    cat("\n")



    #----------------
    # ZIP FILE
    #----------------
if (zip){
    files <- list.files(path=file.path(wd,"output"), recursive=TRUE,full.names = TRUE, include.dirs = TRUE)
    zips <- grep(".zip",files)
    if (length(zips)> 0) {
        files <- files[-zips]
    }

    output <- file.path(wd,"output")
    zip(paste0("BioIndex_results_", paste(as.character(Sys.Date()),format(Sys.time(), "_h%Hm%Ms%OS0"),".zip",sep="")), "output" , root = wd)
    unlink(files , recursive=TRUE)
    unlink(file.path(wd,"output"), force=TRUE, recursive=TRUE)
}


}

