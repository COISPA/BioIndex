#' Main function to perform BioIndex analysis
#'
#' @description
#' BioIndex is an R package designed to support the standardized analysis of MEDITS trawl survey data and the calculation of biological indicators for selected species and population components.
#'
#' @author Walter Zupa \email{zupa@fondazionecoispa.org}
#'
#' @param ta data frame of the TA table in the MEDITS format
#' @param tb data frame of the TB table in the MEDITS format
#' @param tc data frame of the TC table in the MEDITS format
#' @param sspp reference species for the analysis
#' @param rec_threshold cutoff threshold for recruits (reported in mm)
#' @param spaw_threshold cutoff threshold for spawners  (reported in mm)
#' @param haul_threshold minimum number of individuals to be used in estimation of the spatial indicators
#' @param sexes reference sex for the analysis
#' @param depth reference depth range
#' @param GSA reference GSA for the analysis
#' @param country reference country
#' @param map_lim coordinates limits for the maps
#' @param depth_lines depth contours to be plotted in the maps (3 values allowed, e.g c(50,200,800))
#' @param strata data frame of the reference strata for the study area
#' @param stratification_tab data frame of the stratification scheme
#' @param resolution resolution of the depth line
#' @param buffer buffer around the map
#' @param wd path of the working directory
#' @param zip boolean. If TRUE the results are stored in a zip file into the working directory
#' @param save boolean. If TRUE the results are stored in the working directory
#' @param verbose boolean. If TRUE messages are promted in the console
#' @examples
#' \donttest{
#' BioIndex(ta=TA, tb=TB, tc=TC, sspp="MERLMER",rec_threshold=200,
#' spaw_threshold=210,sexes="all", depth=c(10,800), GSA=10, country="all",
#' map_lim=c(13.3,15.2,39.9,41.3),depth_lines=c(50,200,800),
#' strata=BioIndex::strata_scheme, stratification_tab =
#' BioIndex::stratification, resolution=1, buffer=0.1, wd=tempdir(),
#' zip=FALSE, save=TRUE, verbose=TRUE)
#' }
#'
#' @importFrom methods is
#' @importFrom zip zip
#' @return A \code{list} containing the results of the BioIndex workflow, including data frames and plot objects.
#' @export
BioIndex <- function(ta,
                     tb,
                     tc,
                     sspp,
                     rec_threshold,
                     spaw_threshold,
                     haul_threshold = 30,
                     sexes = "all",
                     depth,
                     GSA,
                     country = "all",
                     map_lim,
                     depth_lines = c(10, 200, 800),
                     strata = BioIndex::strata_scheme,
                     stratification_tab = BioIndex::stratification,
                     resolution = NA,
                     buffer = 0.1,
                     wd = NA,
                     zip = TRUE,
                     save=TRUE,
                     verbose = TRUE) {
    if (is.na(wd)) {
        wd <- tempdir()
        if(verbose) message("No directory specified. Results will be saved in the temporary folder: ", wd)
    }

    if (FALSE) {
        # library(BioIndex)
        # wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\Test_BioIndex_package"
        # ta <- read.table(
        #     file.path(wd, "input", "TA GSA18 2017-2020.csv"),
        #     sep = ";",
        #     header = TRUE
        # )
        # tb <- read.table(
        #     file.path(wd, "input", "TB GSA18 2017-2020.csv"),
        #     sep = ";",
        #     header = TRUE
        # )
        # tc <- read.table(
        #     file.path(wd, "input", "TC GSA18 2017-2020.csv"),
        #     sep = ";",
        #     header = TRUE
        # )
        #
        # ta <- read.table(file.path(wd, "input", "TA.csv"),
        #                  sep = ";",
        #                  header = TRUE)
        # tb <- read.table(file.path(wd, "input", "TB.csv"),
        #                  sep = ";",
        #                  header = TRUE)
        # tc <- read.table(file.path(wd, "input", "TC.csv"),
        #                  sep = ";",
        #                  header = TRUE)

        ta <- read.table(file.path(wd,"TA-2020-GRC.csv"),sep=",",header=TRUE)[,-1]
        tb <- read.table(file.path(wd,"TB-2020-GRC.csv"),sep=",",header=TRUE)[,-1]
        tc <- read.table(file.path(wd,"TC-2020-GRC.csv"),sep=",",header=TRUE)[,-1]
        colnames(ta) <- colnames(BioIndex::TA)
        colnames(tb) <- colnames(BioIndex::TB)
        colnames(tc) <- colnames(BioIndex::TC)
        sspp <- "MERLMER"
        rec_threshold = 200
        spaw_threshold = 210
        haul_threshold = 30
        sexes <- "all"
        depth <- c(50, 200)
        GSA <- 20
        country <- "all"
        map_lim <- c(19.166,23.23, 35,39.88)
        depth_lines <- c(50, 100, 200)
        buffer = 0.1
        res = NA
        resolution <- res
        strata = BioIndex::strata_scheme
        stratification_tab = BioIndex::stratification
        save=TRUE
        verbose = TRUE

        BioIndex(
            ta,
            tb,
            tc,
            sspp,
            rec_threshold = rec_threshold,
            spaw_threshold = spaw_threshold,
            sexes = "all",
            depth = depth,
            GSA = GSA,
            country = "all",
            map_lim = map_lim,
            depth_lines = depth_lines,
            strata = strata,
            stratification_tab = stratification_tab,
            resolution = res,
            buffer = buffer,
            wd = wd,
            zip = TRUE,
            save=TRUE,
            verbose = TRUE
        )


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


    if (dir.exists(file.path(wd, "output"))) {
        dir <- list.files(
            file.path(wd, "output"),
            recursive = TRUE,
            full.names = TRUE,
            include.dirs = TRUE
        )
        unlink(dir , recursive = T)
    }



    m <- merge_TATBTC(
        ta,
        tb,
        tc,
        species = sspp,
        country = country,
        wd = wd,
        strata = strata,
        verbose = TRUE
    )

    mTATB <- m[[1]]
    mTATC <- m[[2]]

    ms <- overlayGrid(
        mTATB = mTATB,
        mTATC = mTATC,
        GSA = GSA,
        country = country,
        wd = wd,
        save = save,
        verbose = verbose
    )

    mTATBsp <- ms[[1]]
    mTATCsp <- ms[[2]]

    hauls_position(
        mTATB = mTATB,
        country = country ,
        map_lim = map_lim,
        depth_lines = depth_lines,
        buffer = buffer,
        res = resolution,
        wd = wd,
        save,
        verbose
    )


    #--------------------------------------------------------------
    # Map of abundance and biomass indices by haul
    #--------------------------------------------------------------
    if (verbose) message("\n########################")
    if (verbose) message("Plot of indices by haul")
    if (verbose) message("########################")
    if (verbose) message("")

    source.check <- try(bubble_plot_by_haul_indexes(
        mTATB,
        map_lim,
        depth_lines,
        buffer = 0,
        res = resolution,
        wd = wd,
        save,
        verbose
    )
    ,
    silent = T)

    if (!is(source.check, "try-error"))  {
        if (verbose) message("\nPlot of indices by haul - completed")
    } else {
        if (verbose) message("\nPlot of indices by haul skipped - (Run Error)\n")
    }


    #--------------------------------------------------------------
    # Abundance and biomass indices per GSA in the timeseries
    #--------------------------------------------------------------

    if (verbose) message("\n########################")
    if (verbose) message("Time series of indices")
    if (verbose) message("########################")
    if (verbose) message("")
    index <- indices_ts(
        mTATB,
        GSA = GSA,
        country = country,
        depth_range = depth,
        strata_scheme = strata,
        stratification = stratification_tab,
        wd,
        save
    )
    if (verbose) message("Time series of indices - completed")




    #--------------------------------------------------------------
    # Mean Individual Weights (MIW) per GSA in the timeseries
    #--------------------------------------------------------------

    if (verbose) message("\n########################")
    if (verbose) message("Time series of indices")
    if (verbose) message("########################")
    if (verbose) message("")
    MIW(
        mTATB,
        GSA,
        country,
        depth_range = depth,
        strata_scheme = strata,
        stratification = stratification_tab,
        wd,
        save,
        verbose
    )
    if (verbose) message("Time series of MIW - completed")



    #--------------------------------------------------------------
    # Sex ratio per GSA in the timeseries
    #--------------------------------------------------------------
    if (verbose) message("\n########################")
    if (verbose) message("Sex-ratio time series")
    if (verbose) message("########################")
    if (verbose) message("")

    SR_analysis <- "ok"

    if (SR_analysis == "ok") {
        if (sum(mTATB$NB_OF_FEMALES + mTATB$NB_OF_MALES, na.rm = TRUE) == 0) {
            if (verbose) message("Not enough sex data for sex-ratio estimation")
            if (verbose) message("\nSex-ratio analysis skipped")
        } else {
            sex_ratio(
                mTATB,
                GSA,
                country = country,
                depth_range = depth,
                stratas = strata,
                stratification = stratification_tab,
                wd,
                save
            )
            if (verbose) message("\nSex-ratio analysis - completed")
        }
    } else {
        if (verbose) message("\nSex-ratio analysis skipped")
    }



    # #--------------------------------------------------------------
    # # Abundance indices of spawners in the time series
    # #--------------------------------------------------------------
    # message("\n############################")
    # message("Spawners' abundance indices")
    # message("############################")
    # message("")
    # #------> check the threshold in the file "~/input/maturity_sizes.csv"
    # df_cutoff <- spaw_threshold
    # if (is.na(df_cutoff)) {
    #     message(
    #         "The SPAWNERS' threshold value not provided. Please, define a value in the 'spaw_threshold' parameter."
    #     )
    # }
    # skip_spawners <- FALSE
    # spaw_analysis <- "ok"
    #
    # if (spaw_analysis == "ok") {
    #     source.check <- try(index_spawn(
    #         mTATB,
    #         mTATC,
    #         GSA,
    #         country,
    #         depth_range = depth,
    #         cutoff = spaw_threshold,
    #         stratification = stratification_tab,
    #         wd,
    #         save
    #     )
    #     ,
    #     silent = T)
    #
    #     if (!is(source.check, "try-error"))  {
    #         index_spawn(
    #             mTATB,
    #             mTATC,
    #             GSA,
    #             country,
    #             depth_range = depth,
    #             cutoff = spaw_threshold,
    #             stratification = stratification_tab,
    #             wd,
    #             save
    #         )
    #         if (skip_spawners) {
    #             message("\nSpawners' indices analysis skipped")
    #         } else {
    #             message("\nSpawners' indices analysis - completed")
    #         }
    #
    #     } else {
    #         message("\nSpawners' indices analysis skipped - (Run Error)\n")
    #     }
    #
    # }  else {
    #     message("\nSpawners' indices analysis skipped")
    # }

    #--------------------------------------------------------------
    # Abundance indices of spawners in the time series
    #--------------------------------------------------------------
    if (verbose) message("\n############################")
    if (verbose) message("Spawners' abundance indices")
    if (verbose) message("############################")
    if (verbose) message("")

    df_cutoff <- spaw_threshold
    if (is.na(df_cutoff)) {
        if (verbose) message(
            "The SPAWNERS' threshold value not provided. Please, define a value in the 'spaw_threshold' parameter."
        )
    }

    skip_spawners <- FALSE
    spaw_analysis <- "ok"

    if (spaw_analysis == "ok") {

        source.check <- try(
            index_spawn(
                mTATB,
                mTATC,
                GSA,
                country,
                depth_range = depth,
                cutoff = spaw_threshold,
                stratification = stratification_tab,
                wd,
                save
            ),
            silent = TRUE
        )

        if (!is(source.check, "try-error")) {
            if (skip_spawners) {
                if (verbose) message("\nSpawners' indices analysis skipped")
            } else {
                if (verbose) message("\nSpawners' indices analysis - completed")
            }
        } else {
            if (verbose) message("\nSpawners' indices analysis skipped - (Run Error)\n")
        }

    } else {
        if (verbose) message("\nSpawners' indices analysis skipped")
    }


    # #--------------------------------------------------------------
    # # Abundance indices of recruits in the time series
    # #--------------------------------------------------------------
    # message("\n############################")
    # message("Recruits' abundance indices")
    # message("############################")
    # message("")
    # #------> check the threshold in the file "~/input/maturity_sizes.csv"
    # df_cutoff <- rec_threshold
    # if (is.na(df_cutoff)) {
    #     message(
    #         "The RECRUITS' threshold value not provided. Please, define a value in the 'rec_threshold' parameter."
    #     )
    # }
    # skip_recruits <- FALSE
    # rec_analysis <- "ok"
    #
    # if (rec_analysis == "ok") {
    #     source.check <- try(index_recr(
    #         mTATB,
    #         mTATC,
    #         GSA,
    #         country,
    #         depth_range = depth,
    #         cutoff = rec_threshold,
    #         stratification = stratification_tab,
    #         wd,
    #         save
    #     )
    #     ,
    #     silent = T)
    #
    #     if (!is(source.check, "try-error"))  {
    #         index_recr(
    #             mTATB,
    #             mTATC,
    #             GSA,
    #             country,
    #             depth_range = depth,
    #             cutoff = rec_threshold,
    #             stratification = stratification_tab,
    #             wd,
    #             save
    #         )
    #         if (skip_recruits) {
    #             message("\nRecruits' indices analysis skipped")
    #         } else {
    #             message("\nRecruits' indices analysis - completed")
    #         }
    #
    #     } else {
    #         message("\nRecruits' indices analysis skipped - (Run Error)\n")
    #     }
    #
    # }  else {
    #     message("\nRecruits' indices analysis skipped")
    # }


    #--------------------------------------------------------------
    # Abundance indices of recruits in the time series
    #--------------------------------------------------------------
    if (verbose) message("\n############################")
    if (verbose) message("Recruits' abundance indices")
    if (verbose) message("############################")
    if (verbose) message("")

    df_cutoff <- rec_threshold
    if (is.na(df_cutoff)) {
        if (verbose) message(
            "The RECRUITS' threshold value not provided. Please, define a value in the 'rec_threshold' parameter."
        )
    }

    skip_recruits <- FALSE
    rec_analysis <- "ok"

    if (rec_analysis == "ok") {

        source.check <- try(
            index_recr(
                mTATB,
                mTATC,
                GSA,
                country,
                depth_range = depth,
                cutoff = rec_threshold,
                stratification = stratification_tab,
                wd,
                save
            ),
            silent = TRUE
        )

        if (!is(source.check, "try-error")) {
            if (skip_recruits) {
                if (verbose) message("\nRecruits' indices analysis skipped")
            } else {
                if (verbose) message("\nRecruits' indices analysis - completed")
            }
        } else {
            if (verbose) message("\nRecruits' indices analysis skipped - (Run Error)\n")
        }

    } else {
        if (verbose) message("\nRecruits' indices analysis skipped")
    }


    #--------------------------------------------------------------
    # LFD & L0.95
    #--------------------------------------------------------------
    if (verbose) message("\n############################")
    if (verbose) message("LFD, L0.50 & L0.95")
    if (verbose) message("############################")
    if (verbose) message("")

    LFD_analysis <- "ok"

    if (LFD_analysis == "ok") {
        if (sexes == "all") {
            lfd <- LFD(
                mTATC,
                sex = "all",
                GSA = GSA,
                country = country,
                depth_range = depth,
                strata_scheme = strata,
                stratification = stratification_tab,
                wd,
                save
            )
            Lquant(lfd[[1]], wd, sspp, GSA, save, verbose)
        }

        if (sexes == "M") {
            lfd <- LFD(
                mTATC,
                sex = "M",
                GSA = GSA,
                country = country,
                depth_range = depth,
                strata_scheme = strata,
                stratification = stratification_tab,
                wd,
                save
            )
            Lquant(lfd[[1]], wd, sspp, GSA, save, verbose)
        }
        if (sexes == "F") {
            lfd <- LFD(
                mTATC,
                sex = "F",
                GSA = GSA,
                country = country,
                depth_range = depth,
                strata_scheme = strata,
                stratification = stratification_tab,
                wd,
                save
            )
            Lquant(lfd[[1]], wd, sspp, GSA, save, verbose)
        }
        if (verbose) message("\nLFD & L0.95 analysis - completed")
    } else {
        if (verbose) message("\nLFD & L0.95 analysis skipped")
    }


    #--------------------------------------------------------------
    # Spearman test of trends on short timeseries
    #--------------------------------------------------------------
    if (verbose) message("\n###########################################")
    if (verbose) message("Spearman test of trends on short timeseries")
    if (verbose) message("###########################################")
    if (verbose) message("")
    years <- sort(unique(mTATB$YEAR))
    if (length(years) >= 3) {
        Trend_analysis <- "ok"
        if (Trend_analysis == "ok") {
            spearman(abundance = index[[1]],
                     biomass = index[[2]],
                     years,
                     sspp,
                     wd,
                     save)
            if (verbose) message("\nSpearman test - completed")
        } else {
            if (verbose) message("\nSpearman test skipped")
        }

    } else {
        if (verbose) message("\nSpearman test skipped")
    }


    #--------------------------------------------------------------
    # Abundance indices for statistical squares
    # inverse of CV of abundance indices for statistical squares
    # Biomass indices for statistical squares
    # Mean individual weight for statistical squares
    #--------------------------------------------------------------
    depth_range <- paste(depth, collapse = ",")

    index_on_grid(
        mTATBsp,
        stratum = depth_range,
        wd,
        map_range = map_lim,
        threshold = haul_threshold,
        verbose = TRUE,
        save=TRUE
    )



    #--------------------------------------------------------------
    # Sex ratio for statistical squares
    #--------------------------------------------------------------
    if (verbose) message("\n#############################")
    if (verbose) message("Sex-ratio on GFCM grid")
    if (verbose) message("#############################")
    if (verbose) message("")

    #------> select the minimum number of individuals per haul to be considered in the analysis
    sexratio_grid_skip <- FALSE
    SR_GRID_analysis <- "ok"
    depth_stratum <- paste(depth, collapse = ",")
    if (SR_GRID_analysis == "ok") {
        sex_ratio_on_grid(
            mTATBsp = mTATBsp,
            depth = depth_stratum,
            wd = wd,
            map_range = map_lim,
            threshold = haul_threshold,
            verbose = TRUE,
            save=TRUE
        )

        if (sexratio_grid_skip) {
            if (verbose) message("\nSex-ratio on GFCM grid skipped")
        } else {
            if (verbose) message("\nSex-ratio on GFCM grid - completed")
        }

    } else {
        if (verbose) message("\nSex-ratio on GFCM grid skipped")
    }




    #--------------------------------------------------------------
    # Bubble plots - indices of recruits (abundance)
    # Bubble plots - indices of spawners (abundance)
    #--------------------------------------------------------------
    if (verbose) message("\n################################################")
    if (verbose) message("Bubble plots - indices of recruits and spawners")
    if (verbose) message("################################################")
    if (verbose) message("")
    #------> check the threshold in the file "~/input/maturity_sizes.csv"

    if (any(is.na(c(rec_threshold, spaw_threshold)))) {
        if (verbose) message("Missing threshold value for the selected species.")
    } else {
        bubbleplot_RS_by_hauls(
            mTATC = mTATC,
            map_range = map_lim,
            thresh_rec = rec_threshold,
            thresh_spaw = spaw_threshold,
            depths = depth_lines,
            buffer = buffer,
            res = resolution,
            wd = wd,
            save = save,
            verbose = verbose
        )
        if (verbose) message("\nBubble plots - indices of recruits and spawners - completed")

    }








    if (verbose) message("\n###################")
    if (verbose) message(" Analysis completed")
    if(verbose) message("###################")
    if (verbose) message("")



    #----------------
    # ZIP FILE
    #----------------
    # BioIndex_Report_MERLMER_GSA10_depth_10_800
    if (zip) {
        files <- list.files(
            path = file.path(wd, "output"),
            recursive = TRUE,
            full.names = TRUE,
            include.dirs = TRUE
        )
        zips <- grep(".zip", files)
        if (length(zips) > 0) {
            files <- files[-zips]
        }

        if (length(files) > 0) {
            output <- file.path(wd, "output")
            zip::zip(paste0(
                "BioIndex_results_" ,sspp,"_GSA",GSA,"_Depth",paste(depth, collapse="-"), "m_",
                paste(
                    as.character(Sys.Date()),
                    format(Sys.time(), "_h%Hm%Ms%OS0"),
                    ".zip",
                    sep = ""
                )
            ), "output" , root = wd)
            unlink(files , recursive = TRUE)
            unlink(file.path(wd, "output"),
                   force = TRUE,
                   recursive = TRUE)
            if(verbose) message("ZIP file generated successfully in: ", wd)
        } else {
            if(verbose) message("No files generated to be zipped.")
        }
    }


}
