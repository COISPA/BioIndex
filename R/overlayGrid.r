#' Overlay mTATB and mTATC on GFCM spatial grid
#' @description
#' Overlays survey haul data onto a standard spatial grid (e.g., GFCM grid),
#' harmonizing biological information with a consistent geographic structure
#' for standardized spatial analysis.
#'
#' @param mTATB data frame of the merged TA and TB
#' @param mTATC data frame of the merged TA and TC
#' @param GSA reference GSA for the analysis
#' @param country reference countries for the analysis
#' @param wd working directory used to save results
#' @param save boolean. If TRUE the outputs are saved in the local folder
#' @param verbose boolean. If TRUE messages are prompted in the console
#' @importFrom terra wrap unwrap vect
#' @return A \code{list} or \code{data.frame} containing the merged biological data mapped to the spatial grid.
#' @examples
#' # Use internal data
#' data(TA)
#' data(TB)
#' data(TC)
#' # Create a merged dataset for GSA 10
#' m <- merge_TATBTC(TA[TA$AREA==10,], TB[TB$AREA==10,], TC[TC$AREA==10,], species="MERLMER")
#' mTATB <- m[[1]]
#' mTATC <- m[[2]]
#' overlayGrid(mTATB, mTATC, GSA=10, save=FALSE)
#' @export

overlayGrid <- function(mTATB, mTATC, GSA=NA, country="all", wd=NA, save=TRUE, verbose=FALSE) {
    if (is.na(wd)) { wd <- tempdir() }

    if (FALSE) {
        verbose=TRUE
        save=TRUE
        wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\R_BioIndex_3.3_(in update)"
        ta <- read.table(paste(wd,"input/TA/TA GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
        tb <- read.table(paste(wd,"input/TB/TB GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
        tc <- read.table(paste(wd,"input/TC/TC GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
        country="all"

        # load("D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\- Pacchetto R BioIndex -\\BioIndex\\data\\TA_cols.rda")
        # load("D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\- Pacchetto R BioIndex -\\BioIndex\\data\\TB_cols.rda")
        # load("D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\- Pacchetto R BioIndex -\\BioIndex\\data\\TC_cols.rda")
        species <- "MERLMER"

        m <- merge_TATBTC(ta, tb, tc, species="MERLMER", wd=wd, verbose=TRUE)

        mTATB <- m[[1]]
        mTATC <- m[[2]]

        overlayGrid(mTATB, mTATC, GSA=NA, country="all", wd=wd, save=TRUE, verbose=FALSE)

    }


    if (verbose){
        message("\n########################")
        message("spatial metaDB preparation")
        message("########################")
        message("")
    }

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

    mTATB <- mTATB[mTATB$COUNTRY %in% country_analysis,  ]
    mTATC <- mTATC[mTATC$COUNTRY %in% country_analysis,  ]

    cgpmgrid <- terra::unwrap(BioIndex::cgpmgrid)

    metaDB <- mTATB
    if (any(is.na(GSA))) {
        GSA <- unique(metaDB$GSA)[1]
    }
    metaDB <- metaDB[metaDB$GSA %in% GSA, ]
    GENERE <- as.character(unique(metaDB$GENUS)[unique(metaDB$GENUS) != -1])
    SPECIE <- as.character(unique(metaDB$SPECIES)[unique(metaDB$SPECIES) != -1])
    sspp <- paste(GENERE,SPECIE, sep="")

    # Overlay metaDB with the grid
    dft <- data.frame(metaDB[, c("MEAN_LONGITUDE_DEC","MEAN_LATITUDE_DEC")])
    vt <- terra::vect(dft, geom=c("MEAN_LONGITUDE_DEC", "MEAN_LATITUDE_DEC"), crs="+proj=longlat", keepgeom=FALSE)
    overlay <- terra::extract(x=cgpmgrid, y=vt)
    dupl <- which(duplicated(overlay$id.y))
    if (length(dupl)>0){
        overlay <- overlay[-(which(duplicated(overlay$id.y))), ]
    }

    metaDBnew_georef <- data.frame(metaDB, x_center= overlay$gfcm.cen_5,y_center=overlay$gfcm.cen_6 ,cgpmgrid_id=overlay$GFCM_ID)

    centroidi <- terra::unwrap(BioIndex::centroidi)

    mTATB_geo <- metaDBnew_georef
    if(verbose){
        message("Catch metaDB successfully created.")
    }
    if (save){

        write.table(metaDBnew_georef, paste(wd, "/output/",sspp," - allGSAs_metaDB_catch in GRID.csv",sep=""), sep=";", row.names=F)
        if(verbose){
            message(paste("Catch metaDB saved in the output folder"))
            message("")
        }
    }
    if (any(is.na(metaDBnew_georef$cgpmgrid_id))) {
        if (verbose) {
            message("One or more hauls are out of the reference grid area. Data from these hauls should be excluded by statial analysis.")
            message("")
        }
    }

    # Read metaDB file and transform into spatial object
    metaDB <- mTATC
    metaDB <- metaDB[metaDB$GSA %in% GSA, ]
    overlay <- terra::extract(cgpmgrid, data.frame(metaDB$MEAN_LONGITUDE_DEC, metaDB$MEAN_LATITUDE_DEC))

    dupl <- NULL
    dupl <- which(duplicated(overlay$id.y))
    if (length(dupl)>0){
        overlay <- overlay[-(which(duplicated(overlay$id.y))), ]
    }
    metaDBnew_georef <- data.frame(metaDB, x_center= overlay$gfcm.cen_5,y_center=overlay$gfcm.cen_6 ,cgpmgrid_id=overlay$GFCM_ID)

    mTATC_geo <- metaDBnew_georef
    if(verbose){
        message("Biological metaDB successfully created.")
    }
    if(save){

        write.table(metaDBnew_georef, paste(wd, "/output/",sspp," - allGSAs_metaDB_biological in GRID.csv",sep=""), sep=";", row.names=F)
        if(verbose){
            message(paste("Biological metaDB saved in the output folder"))
            message("")
        }
    }
    if (any(is.na(metaDBnew_georef$cgpmgrid_id))) {
        if (verbose){
            message("One or more hauls are out of the reference grid area. Data from these hauls should be excluded by statial analysis.")
            message("")
        }
    }
    return(list(mTATB_geo,mTATC_geo))
}