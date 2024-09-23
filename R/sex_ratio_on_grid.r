#' Plot sex ratio spatial distribution
#'
#' @param mTATBsp spatial mTATB
#' @param depth reference depth range
#' @param wd working directory
#' @param map_range range of coordinates for the map
#' @param threshold minimum number of individuals per haul
#' @param save boolean. If TRUE the outputs are saved in the local folder
#' @param verbose boolean. If TRUE messages are prompted in the console
#' @importFrom terra extract unwrap
#' @importFrom stats quantile
#' @export

sex_ratio_on_grid <- function(mTATBsp, depth, wd, map_range,threshold=30,verbose=FALSE, save=FALSE) {

    if (FALSE) {
        verbose=TRUE
        save=TRUE
        map_range <- c(15,21,39.5,42.5)
        wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\R_BioIndex_3.3_(in update)"
        ta <- read.table(paste(wd,"input/TA/TA GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
        tb <- read.table(paste(wd,"input/TB/TB GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
        tc <- read.table(paste(wd,"input/TC/TC GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
        country="all"
        species <- "MERLMER"

        m <- merge_TATBTC(ta, tb, tc, species="MERLMER", country="all", wd=wd, verbose=TRUE)

        mTATB <- m[[1]]
        mTATC <- m[[2]]

        mTATBsp <- overlayGrid(mTATB, mTATC, GSA=NA, wd=wd, save=TRUE, verbose=FALSE)[[1]]
        threshold=30
        depth <- "10,800"
        sex_ratio_on_grid (mTATBsp, depth, wd, map_range,threshold=30,verbose=TRUE, save=TRUE)
    }

    level <- category <- long <- lat <- group <- NULL

    depthstr <- depth
    depth_range <- as.numeric(str_split(depth,",")[[1]])

    if (depthstr == "5,45") {
        cgpmgrid <- terra::unwrap(BioIndex::stratum_0_45)
    }

    if (depthstr == "5,35") {
        cgpmgrid <- terra::unwrap(BioIndex::stratum_0_35)
    }

    if (depthstr == "10,125") {
        cgpmgrid <- terra::unwrap(BioIndex::stratum_0_125)
    }

    if (depthstr == "5,125") {
        cgpmgrid <- terra::unwrap(BioIndex::stratum_0_125)
    }

    if (depthstr == "0,125") {
        cgpmgrid <- terra::unwrap(BioIndex::stratum_0_125)
    }

    if (depthstr == "10,200") {
        cgpmgrid <- terra::unwrap(BioIndex::stratum_0_200)
    }

    if (depthstr == "10,800") {
        cgpmgrid <- terra::unwrap(BioIndex::stratum_0_800)
    }

    if (depthstr == "200,800") {
        cgpmgrid <- terra::unwrap(BioIndex::stratum_200_800)
    }
    cgpmgrid_bkp <- cgpmgrid
    continent <- terra::unwrap(BioIndex::continent)



    metaDBnew <- mTATBsp
    years <- sort(unique(metaDBnew$YEAR))
    GSA <- unique(metaDBnew$GSA)[1]
    GENERE <- as.character(unique(metaDBnew$GENUS)[unique(metaDBnew$GENUS) != -1])
    SPECIE <- as.character(unique(metaDBnew$SPECIES)[unique(metaDBnew$SPECIES) != -1])
    sspp <- paste(GENERE,SPECIE, sep="")

    if (length(years) < 10) {
        n_years <- length(years)
        firstyear <- years[length(years)-(n_years-1)]
        last_10year <- years[years >= firstyear]
    } else {
        firstyear <- years[length(years)-9]
        last_10year <- years[years >= firstyear]
    }

    metaDBnew <- metaDBnew[metaDBnew$YEAR %in% last_10year & metaDBnew$TOTAL_NUMBER_IN_THE_HAUL > threshold, ]  #
    metaDBnew <- metaDBnew[metaDBnew$MEAN_DEPTH > depth_range[1] & metaDBnew$MEAN_DEPTH <= depth_range[2], ]  #


    if (nrow(metaDBnew)==0){
        if (verbose){
        cat("Not enough data for the sex ratio estimation.\nTry again reducing the threshold value.\n")
            }
        sexratio_grid_skip <- TRUE
    } else {

        if ( all(metaDBnew$NB_OF_FEMALES == 0) & all(metaDBnew$NB_OF_MALES == 0) ) {
            if(verbose){
            cat("No sex information available in TC for the selected species. Sex-ratio estimation is not possible")
                }
            sexratio_grid_skip <- TRUE
        } else {

            grid_id <- unique(metaDBnew$cgpmgrid_id)
            sex_ratio <- data.frame(matrix(NA, nrow=length(grid_id), ncol=4))
            colnames(sex_ratio) <- c("cgpmgrid_id", "ratio", "sd", "CV")
            sex_ratio[,1] <-grid_id

            i=3
            for (i in 1:length(grid_id)){
                dd <- metaDBnew[metaDBnew$cgpmgrid_id == grid_id[i], ]
                dd$nkm2_F <- dd$NB_OF_FEMALES/dd$SWEPT_AREA
                dd$nkm2_M <- dd$NB_OF_MALES/dd$SWEPT_AREA
                dd$nkm_FM <- (dd$NB_OF_MALES + dd$NB_OF_FEMALES)/dd$SWEPT_AREA
                dd$sr_haul<- dd$nkm2_F / dd$nkm_FM

                nkm_F <- sum(dd$nkm2_F)
                nkm_M <- sum(dd$nkm2_M)
                nkm_FM <- sum(dd$nkm_FM)
                sex_ratio[i,2] <- nkm_F/nkm_FM
                sex_ratio[i,3] <- sqrt(sum((sex_ratio[i,2] - dd$sr_haul)^2)/length(dd$sr_haul))
                sex_ratio[i,4] <- sex_ratio[i,3]/sex_ratio[i,2]
            }

            if(save){
            write.table(sex_ratio, paste(wd, "/output/",sspp," - GFCM SEX RATIO.csv", sep=""), sep=";", row.names=FALSE)
            }

            names(map_range) <- c("xmin","xmax","ymin","ymax")
            x1 <- min(map_range[1])
            x2 <- max(map_range[2])
            y1 <- min(map_range[3])
            y2 <- max(map_range[4])

            #----------------------------------------
            # Saving maps of the sex-ratio by grid
            #----------------------------------------
            cat_lim <- quantile(sex_ratio$ratio[!is.na(sex_ratio$ratio)], probs = seq(0, 1,0.20) )
            cate_label <-  round(cat_lim,3)

            cate_label_min <- cate_label
            cate_label_max <- cate_label[2:length(cate_label)]

            cate_label <- paste("(", cate_label_min, ", ", cate_label_max, "]", sep="")
            cate_label[1] <- paste("[", cate_label_min[1], ", ", cate_label_max[1], "]", sep="")

            cate_label <- cate_label[1: (length(cate_label)-1)]

            palette_cate <- heat.colors(length(cate_label))
            palette_cate <- palette_cate[length(palette_cate):1]

            # continent <- unwrap(BioIndex::continent)
            cgpmgrid <- cgpmgrid_bkp
            cgpmgrid$SR <- NA
            cgpmgrid$color <- "#FFFFFF"

            i=1
            for (i in 1:length(sex_ratio[,1])){
                cgpmgrid$SR[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] <- sex_ratio[i,"ratio"]

                if (all(is.na(cgpmgrid$SR[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]]) )){
                    cgpmgrid$color[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] <- "#FFFFFF"
                } else {

                    if(all(cgpmgrid$SR[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] <= cat_lim[2])  & all(cgpmgrid$SR[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] >= cat_lim[1])){
                        cgpmgrid$color[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] <- palette_cate[1]
                    }

                    if(all(cgpmgrid$SR[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] <= cat_lim[3])  & all(cgpmgrid$SR[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] > cat_lim[2])){
                        cgpmgrid$color[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] <- palette_cate[2]
                    }

                    if(all(cgpmgrid$SR[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] <= cat_lim[4])  & all(cgpmgrid$SR[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] > cat_lim[3])){
                        cgpmgrid$color[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] <- palette_cate[3]
                    }

                    if(all(cgpmgrid$SR[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] <= cat_lim[5])  & all(cgpmgrid$SR[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] > cat_lim[4])){
                        cgpmgrid$color[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] <- palette_cate[4]
                    }

                    if(all(cgpmgrid$SR[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] <= cat_lim[6])  & all(cgpmgrid$SR[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] > cat_lim[5])){
                        cgpmgrid$color[cgpmgrid$GFCM_ID==sex_ratio[i,"cgpmgrid_id"]] <- palette_cate[5]
                    }
                }
            }


            lx =map_range["xmax"] - map_range["xmin"]
            ly =map_range["ymax"] - map_range["ymin"]
            ratio <- ly/lx
            img_width <- 30
            img_height <- img_width*ratio


            # par( mar=c(4, 5, 4, 2)) #c(bottom, left, top, right)
            # plot(1,1,type="n",xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), xlab=expression(paste("Longitude (",degree,"E)")), ylab=expression(paste("Latitude (",degree,"N)")), main=paste("sex ratio -", sspp, " GSA", GSA),asp=1 )
            # plot(cgpmgrid, col=cgpmgrid$color, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")                                                                                 #data = centroidi_coords,
            # if (depth == "200,800" | depth == 3){
            #     plot(mask,col="white", xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")
            # }
            # plot(continent, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), border="grey", col="light grey", add=T)
            # legend("topleft", cate_label, col=palette_cate, pch = 15, pt.cex=c(2), title = "sex ratio\n(classes based on percentiles)", bty = "n")
            # box()
            # psr <- recordPlot()

            tab_color <- data.frame(category=cate_label,color=palette_cate,level=seq(1,length(cate_label),1))
            cgpmgrid <- terra::merge(cgpmgrid,tab_color,by.x="color",by.y="color")
            world <- map_data("world")
            xl <- c(x1,x2)
            yl <- c(y1,y2)
            x_breaks <- c(round(xl[1], 0), round(xl[1], 0) + round((xl[2] - xl[1]) / 2, 0), round(xl[1], 0) + 2 * round((xl[2] - xl[1]) / 2, 0))
            y_breaks <- c(round(yl[1], 0), round(yl[1], 0) + round((yl[2] - yl[1]) / 2, 0), round(yl[1], 0) + 2 * round((yl[2] - yl[1]) / 2, 0))
            buff=0

            cols <- c("#FFFF80",
                      "#FFFF00",
                      "#FFAA00",
                      "#FF5500",
                      "#FF0000")

            names(cols) <- tab_color$category
            cgpmgrid$level <- is.factor(cgpmgrid$level)
            cgpmgrid$level <- factor(cgpmgrid$level, levels=seq(1,nrow(tab_color),1))
            cgpmgrid <- cgpmgrid %>% arrange(level)

            p1 <- ggplot(cgpmgrid)+
                geom_spatvector(aes(fill = category))+
                scale_fill_manual(values = cols, breaks=names(cols))+
                geom_polygon(data=world, aes(long,lat,group=group), fill="dark grey", colour="darkgrey")+
                coord_sf(xlim = xl, ylim = yl)+
                labs(fill="Sex ratio")+
                xlab("longitude")+
                ylab("latitude")



            print(p1)

            if (save){
                jpeg(filename=paste(wd, "/output/",sspp," - GFCM GRID Sex Ratio.jpg", sep = ""), width=25, height=25, bg="white", units="cm",res=200)
                print(p1)
                dev.off()
            }

            return(sex_ratio)
        }
    }

}
