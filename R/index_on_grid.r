#' Generating maps of indexes
#'
#' @param mTATBsp spatial mTATB
#' @param stratum reference stratum range (allowed values: "10,200","10,800","200,800","5,35","5,45")
#' @param wd working directory
#' @param map_range range of coordinates for the map
#' @param threshold minimum number of individuals per haul
#' @param verbose boolean. If TRUE messages are promted in the console
#' @param save boolean. If TRUE the results are stored in the working directory
#' @export
#' @importFrom ggplot2 ggplot scale_fill_manual geom_polygon coord_sf labs xlab ylab map_data
#' @importFrom stringr str_split
#' @importFrom tidyterra geom_spatvector
#' @importFrom dplyr arrange
#' @importFrom stats quantile

index_on_grid <- function(mTATBsp, stratum, wd, map_range, threshold = 30, verbose = FALSE, save = FALSE) {
  if (FALSE) {
    # library(stringr)
    # library(ggplot2)
    # library(tidyterra)
    verbose <- TRUE
    save <- TRUE
    map_range <- c(15, 21, 39.5, 42.5)
    wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\Test_BioIndex_package"
    ta <- read.table(paste(wd, "input/TA GSA18 2017-2020.csv", sep = "/"), sep = ";", header = T)
    tb <- read.table(paste(wd, "input/TB GSA18 2017-2020.csv", sep = "/"), sep = ";", header = T)
    tc <- read.table(paste(wd, "input/TC GSA18 2017-2020.csv", sep = "/"), sep = ";", header = T)
    country <- "all"
    species <- "MERLMER"

    m <- merge_TATBTC(ta, tb, tc, species = "MERLMER", country = "all", wd = wd, verbose = TRUE)

    mTATB <- m[[1]]
    mTATC <- m[[2]]

    mTATBsp <- overlayGrid(mTATB, mTATC, GSA = NA, wd = wd, save = TRUE, verbose = FALSE)[[1]]

    stratum <- "10,800"
    index_on_grid(mTATBsp, stratum, wd, map_range, threshold = 30, verbose = FALSE, save = TRUE)
  }

  level <- category <- long <- lat <- group <- NULL

  # threshold_to_ind_weight_calc <- 30

  depthstr <- stratum
  depth_range <- as.numeric(str_split(stratum,",")[[1]])

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

  metaDB <- mTATBsp
  metaDB <- metaDB[metaDB$MEAN_DEPTH > depth_range[1] & metaDB$MEAN_DEPTH <= depth_range[2], ]

  # }
  GENERE <- as.character(unique(metaDB$GENUS)[unique(metaDB$GENUS) != -1])
  SPECIE <- as.character(unique(metaDB$SPECIES)[unique(metaDB$SPECIES) != -1])
  sspp <- paste(GENERE, SPECIE, sep = "")

  years <- sort(unique(metaDB[, "YEAR"]))

  if (length(years) < 10) {
    n_years <- length(years)
    firstyear <- years[length(years) - (n_years - 1)]
    last_10year <- years[years >= firstyear]
  } else {
    firstyear <- years[length(years) - 9]
    last_10year <- years[years >= firstyear]
  }

  metaDB <- metaDB[metaDB$YEAR %in% last_10year, ]

  colnames(metaDB)[colnames(metaDB) == "cgpmgrid_id"] <- "CGPMgrid"
  colnames(metaDB)[colnames(metaDB) == "kg_km2"] <- "B_km2"
  colnames(metaDB)[colnames(metaDB) == "kg_h"] <- "B_h"

  # Create Mean Fish Weight by haul
  metaDB$MEAN_WEIGHT <- round(metaDB$TOTAL_WEIGHT_IN_THE_HAUL / metaDB$TOTAL_NUMBER_IN_THE_HAUL, 3) / 1000

  #----------------------------------------
  # CALCULATION OF THE ABUNDANCE and saving the tables
  #----------------------------------------
  nhauls.grid <- aggregate(metaDB$N_km2, by = list(metaDB$CGPMgrid, metaDB$x_center, metaDB$y_center), length) # [,4]
  colnames(nhauls.grid) <- c("CGPMgrid", "x_center", "y_center", "nhauls")

  cgpmgridlevel <- unique(metaDB$CGPMgrid)
  ntab <- data.frame(CGPMgrid = cgpmgridlevel, GSA = unique(metaDB$GSA))
  ntab$meanNkm2 <- NA
  i <- 1
  for (i in 1:length(cgpmgridlevel)) {
    ddd <- metaDB[metaDB$CGPMgrid == cgpmgridlevel[i], ]
    ntab[i, 3] <- mean(ddd[!is.na(ddd$N_km2), "N_km2"])
    ssd <- (ddd$N_km2 - ntab[i, 3])^2 # sum of squared differences
    variance <- sum(ssd) / (length(ddd$N_km2) - 1) # variance
    ntab$sdNkm2[i] <- sqrt(variance)
    ntab$cvNkm2[i] <- sqrt(variance) / ntab[i, 3]
    ntab$inverse_cvNkm2[i] <- 1 / ntab$cvNkm2[i]
    ntab$nhauls_in_square[i] <- length(ddd[, "HAUL_NUMBER"])
    ntab$nhauls_positive[i] <- length(ddd[ddd$N_km2 > 0, "HAUL_NUMBER"])
  }

  summary.grid <- merge(ntab, nhauls.grid[, 1:3])
  colnames(summary.grid) <- c("cgpmgridlevel", "GSA", "meanNkm2", "sdNkm2", "cvNkm2", "inverse_cvNkm2", "nhauls", "positive_hauls", "lon", "lat")

  abundance_grid <- summary.grid
  abundance_grid$sspp <- sspp
  if (save) {
    write.table(abundance_grid, paste(wd, "/output/", sspp, " - GFCM GRID ABUNDANCE.csv", sep = ""), sep = ";", row.names = FALSE)

    if (verbose) {
      cat("Abundance indices for statistical squares correctly estimated \ninverse of CV of abundance indices for statistical squares correctly estimated \n")
      cat(paste("file of abundance indices for statistical squares saved in the following folder: '", wd, "/output/", sspp, " - GFCM GRID ABUNDANCE.csv \n", sep = ""))
    }
  }
  #----------------------------------------
  # CALCULATION OF THE BIOMASS and saving the tables
  #----------------------------------------
  nhauls.grid <- aggregate(metaDB$B_km2, by = list(metaDB$CGPMgrid, metaDB$x_center, metaDB$y_center), length) # [,4]
  colnames(nhauls.grid) <- c("CGPMgrid", "x_center", "y_center", "nhauls")

  cgpmgridlevel <- unique(metaDB$CGPMgrid)
  Btab <- data.frame(CGPMgrid = cgpmgridlevel, GSA = unique(metaDB$GSA))
  Btab$meanBkm2 <- NA
  i <- 1
  for (i in 1:length(cgpmgridlevel)) {
    ddd <- metaDB[metaDB$CGPMgrid == cgpmgridlevel[i], ]
    Btab[i, 3] <- mean(ddd[!is.na(ddd$B_km2), "B_km2"])
    ssd <- (ddd$B_km2 - Btab[i, 3])^2 # sum of squared differences
    variance <- sum(ssd) / (length(ddd$B_km2) - 1) # variance
    Btab$sdBkm2[i] <- sqrt(variance)
    Btab$cvBkm2[i] <- sqrt(variance) / Btab[i, 3]
    Btab$inverse_cvBkm2[i] <- 1 / Btab$cvBkm2[i]
    Btab$nhauls_in_square[i] <- length(ddd[, "HAUL_NUMBER"])
    Btab$nhauls_positive[i] <- length(ddd[ddd$B_km2 > 0, "HAUL_NUMBER"])
  }

  summary.grid <- merge(Btab, nhauls.grid[, 1:3])
  colnames(summary.grid) <- c("cgpmgridlevel", "GSA", "meanBkm2", "sdBkm2", "cvBkm2", "inverse_cvBkm2", "nhauls", "positive_hauls", "lon", "lat")

  biomass_grid <- summary.grid
  biomass_grid$sspp <- sspp
  if (save) {
    write.table(biomass_grid, paste(wd, "/output/", sspp, " - GFCM GRID BIOMASS.csv", sep = ""), sep = ";", row.names = FALSE)
    if (verbose) {
      cat("Biomass indices for statistical squares correctly estimated \n")
      cat(paste("file of Biomass indices for statistical squares saved in the following folder: '", wd, "/output/", sspp, " - GFCM GRID BIOMASS.csv \n", sep = ""))
    }
  }
  #----------------------------------------
  # Calculating mean weight
  #----------------------------------------

  nhauls.grid <- aggregate(metaDB$N_km2, by = list(metaDB$CGPMgrid, metaDB$x_center, metaDB$y_center), length) # [,4]
  colnames(nhauls.grid) <- c("CGPMgrid", "x_center", "y_center", "nhauls")

  cgpmgridlevel <- unique(metaDB$CGPMgrid)
  MIWtab <- data.frame(CGPMgrid = cgpmgridlevel, GSA = unique(metaDB$GSA))
  MIWtab$MIW <- NA
  i <- 1
  for (i in 1:length(cgpmgridlevel)) {
    ddd <- metaDB[metaDB$CGPMgrid == cgpmgridlevel[i], ]
    MIWtab[i, 3] <- sum(ddd[!is.na(ddd$MEAN_WEIGHT) & ddd$MEAN_WEIGHT > 0, "MEAN_WEIGHT"]) / length(ddd[!is.na(ddd$MEAN_WEIGHT) & ddd$MEAN_WEIGHT > 0, "MEAN_WEIGHT"])
    ssd <- (ddd[!is.na(ddd$MEAN_WEIGHT) & ddd$MEAN_WEIGHT > 0, "MEAN_WEIGHT"] - MIWtab[i, 3])^2 # sum of squared differences
    variance <- sum(ssd) / (length(ddd[!is.na(ddd$MEAN_WEIGHT) & ddd$MEAN_WEIGHT > 0, "MEAN_WEIGHT"]) - 1) # variance
    MIWtab$sdNkm2[i] <- sqrt(variance)
    MIWtab$cvNkm2[i] <- sqrt(variance) / MIWtab[i, 3]
    MIWtab$inverse_cvNkm2[i] <- 1 / MIWtab$cvNkm2[i]
    MIWtab$nhauls_in_square[i] <- length(ddd[, "HAUL_NUMBER"])
    MIWtab$nhauls_positive[i] <- length(ddd[!is.na(ddd$N_km2) & ddd$N_km2 > 0, "HAUL_NUMBER"])
  }

  summary.grid <- merge(MIWtab, nhauls.grid[, 1:3])
  colnames(summary.grid) <- c("cgpmgridlevel", "GSA", "MIW", "sdMIW", "cvMIW", "inverse_cvMIW", "nhauls", "positive_hauls", "lon", "lat")

  meanWEIGHT_grid <- summary.grid
  meanWEIGHT_grid$sspp <- sspp
  if (save) {
    write.table(meanWEIGHT_grid, paste(wd, "/output/", sspp, " - GFCM GRID MIW.csv", sep = ""), sep = ";", row.names = FALSE)
    if (verbose) {
      cat("MIW for statistical squares correctly estimated \ninverse of CV of MIW for statistical squares correctly estimated \n")
      cat(paste("file of MIW for statistical squares saved in the following folder:\n '", wd, "/output/", sspp, " - GFCM GRID MIW.csv \n", sep = ""))
    }
  }

  #----------------------------------------
  # Plotting indices
  #----------------------------------------

  # plotGFCM(
  #   abundance_grid = abundance_grid,
  #   biomass_grid = biomass_grid,
  #   meanWEIGHT_grid = meanWEIGHT_grid,
  #   map_range = map_range,
  #   stratum = stratum
  # )



  sspp <- unique(abundance_grid$sspp)[1]
  names(map_range) <- c("xmin","xmax","ymin","ymax")
  summary.grid <- abundance_grid
  # metaDB <- mTATB

  x1 <- min(map_range[1])
  x2 <- max(map_range[2])
  y1 <- min(map_range[3])
  y2 <- max(map_range[4])

  #----------------------------------------
  # Saving maps of the ABUNDANCE index
  #----------------------------------------                                                                          # [summary.grid$meanNkm2>5]
  cat_lim <- as.numeric(quantile(summary.grid[!is.na(summary.grid$meanNkm2), "meanNkm2"], probs = seq(0, 1,0.20) ) )
  cate_label <-   round(cat_lim,1)
  cate_label_min <- cate_label
  cate_label_max <- cate_label[2:length(cate_label)]


  cate_label <- paste("(", cate_label_min, ", ", cate_label_max, "]", sep="")
  cate_label[1] <- paste("[", cate_label_min[1], ", ", cate_label_max[1], "]", sep="")
  cate_label <- cate_label[1: (length(cate_label)-1)]

  palette_cate <- heat.colors(length(cate_label))
  palette_cate <- palette_cate[length(palette_cate):1]

  # cgpmgrid <- terra::unwrap(cgpmgrid)
  # continent <- terra::unwrap(continent)
  cgpmgrid$abundance <- NA
  cgpmgrid$color <- "#FFFFFF"
  i=2
  for (i in 1:length(summary.grid[,1])){
    cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- summary.grid[i,"meanNkm2"]
    if(all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[2])  & all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] >= cat_lim[1])){
      cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[1]
    }

    if(all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[3] ) & all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[2])){
      cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[2]
    }

    if(all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[4])  & all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[3])){
      cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[3]
    }

    if(all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[5])  & all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[4])){
      cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[4]
    }

    if(all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[6])  & all(cgpmgrid$abundance[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[5])){
      cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[5]
    }
  }


  lx =map_range["xmax"] - map_range["xmin"]
  ly =map_range["ymax"] - map_range["ymin"]
  ratio <- ly/lx
  img_width <- 30
  img_height <- img_width*ratio

  # par( mar=c(4, 5, 4, 2)) #c(bottom, left, top, right)
  # plot(1,1,type="n",xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]) ,  xlab=expression(paste("Longitude (",degree,"E)")), ylab=expression(paste("Latitude (",degree,"N)")), main=paste(paste("Abundance  (n/km^2)", sep=" "), sspp), asp=1)
  # plot(cgpmgrid,col=cgpmgrid$color, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=TRUE, border="light grey")                                                                                 #data = centroidi_coords,
  # if (stratum == "200,800" | stratum == 3){
  #   plot(mask,col="white", xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")
  # }
  # plot(continent, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), border="grey", col="light grey", add=T)
  # legend("topleft", cate_label, col=palette_cate, pch = 15, pt.cex=c(2), title = "Abundance (n/km^2)\n(classes based on percentiles)", bty = "n")
  # # legend(map_range[1,"xmin"],map_range[1,"ymax"], cate_label, col=palette_cate, pch = 15, pt.cex=c(2), title = "Abundance (n/km^2)\n(classes based on percentiles)", bty = "n")
  # box()
  # p1 <- recordPlot()
  #
  # jpeg(file=paste(wd, "/output/",sspp," - GFCM GRID ABUNDANCE.jpg", sep = ""), width=25, height=25, bg="white", units="cm",res=200)
  # replayPlot(p1)
  # dev.off()

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
 # cgpmgrid$level <- as.character( cgpmgrid$level)
 cgpmgrid$level <- is.factor(cgpmgrid$level)
 cgpmgrid$level <- factor(cgpmgrid$level, levels=seq(1,nrow(tab_color),1))
 cgpmgrid <- cgpmgrid %>% arrange(level)

 p1 <- ggplot(cgpmgrid)+
   geom_spatvector(aes(fill = category))+
   scale_fill_manual(values = cols, breaks=names(cols))+
   geom_polygon(data=world, aes(long,lat,group=group), fill="dark grey", colour="darkgrey")+
   coord_sf(xlim = xl, ylim = yl)+
   labs(fill="Abundance\n (n/km^2)")+
   xlab("longitude")+
   ylab("latitude")

 print(p1)

 if (save){
   jpeg(filename=paste(wd, "/output/",sspp," - GFCM GRID ABUNDANCE.jpg", sep = ""), width=25, height=25, bg="white", units="cm",res=200)
   print(p1)
   dev.off()
 }
  #----------------------------------------
  # Saving maps of the inverse of CV of ABUNDANCE index
  #----------------------------------------

  cgpmgrid <- cgpmgrid_bkp
  cat_lim <- quantile(summary.grid$inverse_cvNkm2[!is.na(summary.grid$inverse_cvNkm2)], probs = seq(0, 1,0.20) )
  cate_label <-  round(cat_lim,3)

  cate_label_min <- cate_label
  cate_label_max <- cate_label[2:length(cate_label)]


  cate_label <- paste("(", cate_label_min, ", ", cate_label_max, "]", sep="")
  cate_label[1] <- paste("[", cate_label_min[1], ", ", cate_label_max[1], "]", sep="")

  #cate_label[length(cate_label)] <- paste(">", cate_label_max[length(cate_label_max)] )
  cate_label <- cate_label[1: (length(cate_label)-1)]

  palette_cate <- heat.colors(length(cate_label))
  palette_cate <- palette_cate[length(palette_cate):1]

  cgpmgrid$invCV <- NA
  cgpmgrid$color <- "#FFFFFF"

  i=3
  for (i in 1:length(summary.grid[,1])){
    cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- summary.grid[i,"inverse_cvNkm2"]

    if(all(is.na(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]])) ){
      cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- "#FFFFFF"
    } else {

      if(all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[2])  & all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] >= cat_lim[1])){
        cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[1]
      }

      if(all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[3])  & all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[2])){
        cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[2]
      }

      if(all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[4])  & all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[3])){
        cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[3]
      }

      if(all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[5])  & all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[4])){
        cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[4]
      }

      if(all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[6])  & all(cgpmgrid$invCV[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[5])){
        cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[5]
      }
    }
  }

  lx =map_range["xmax"] - map_range["xmin"]
  ly =map_range["ymax"] - map_range["ymin"]
  ratio <- ly/lx
  img_width <- 30
  img_height <- img_width*ratio

  # par( mar=c(4, 5, 4, 2)) #c(bottom, left, top, right)
  # plot(1,1,type="n",xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), xlab=expression(paste("Longitude (",degree,"E)")), ylab=expression(paste("Latitude (",degree,"N)")), main=paste("Inverse of Abundance CV -", sspp) , asp=1)
  # plot(cgpmgrid,col=cgpmgrid$color, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")                                                                                 #data = centroidi_coords,
  # if (stratum == "200,800" | stratum == 3){
  #   plot(mask,col="white", xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")
  # }
  # plot(continent, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), border="grey", col="light grey", add=T)
  # legend("topleft", cate_label, col=palette_cate, pch = 15, pt.cex=c(2), title = "Inverse of Abundance CV\n(classes based on percentiles)", bty = "n")
  # box()
  # p2 <- recordPlot()

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
  # cgpmgrid$level <- as.character( cgpmgrid$level)
  cgpmgrid$level <- is.factor(cgpmgrid$level)
  cgpmgrid$level <- factor(cgpmgrid$level, levels=seq(1,nrow(tab_color),1))
  cgpmgrid <- cgpmgrid %>% arrange(level)

  p2 <- ggplot(cgpmgrid)+
    geom_spatvector(aes(fill = category))+
    scale_fill_manual(values = cols, breaks=names(cols))+
    geom_polygon(data=world, aes(long,lat,group=group), fill="dark grey", colour="darkgrey")+
    coord_sf(xlim = xl, ylim = yl)+
    labs(fill="Inverse of\n Abundance CV")+
    xlab("longitude")+
    ylab("latitude")

  print(p2)

  if(save){
    jpeg(filename=paste(wd, "/output/",sspp," - GFCM GRID ABUNDANCE Inverse CV.jpg", sep = ""), width=25, height=25, bg="white", units="cm",res=200)
    print(p2)
    dev.off()
  }


  ###################
  ###################
  ###             ###
  ###   BIOMASS   ###
  ###             ###
  ###################
  ###################

  cgpmgrid <- cgpmgrid_bkp
  summary.grid <- biomass_grid
  # metaDB <- mTATB

  #----------------------------------------
  # Saving maps of the BIOMASS index
  #----------------------------------------

  cat_lim <- quantile(summary.grid[!is.na(summary.grid$meanBkm2), "meanBkm2"], probs = seq(0, 1,0.20) )
  cate_label <-  round(cat_lim,3)

  cate_label_min <- cate_label
  # cate_label_min[1] <- 0
  cate_label_max <- cate_label[2:length(cate_label)]


  cate_label <- paste("(", cate_label_min, ", ", cate_label_max, "]", sep="")
  cate_label[1] <- paste("[", cate_label_min[1], ", ", cate_label_max[1], "]", sep="")

  cate_label <- cate_label[1: (length(cate_label)-1)]
  # cate_label[length(cate_label)] <- paste(">", cate_label_max[length(cate_label_max)] )


  palette_cate <- heat.colors(length(cate_label))
  palette_cate <- palette_cate[length(palette_cate):1]


  cgpmgrid$meanBkm2 <- NA
  cgpmgrid$color <- "#FFFFFF"

  i=3
  for (i in 1:length(summary.grid[,1])){
    cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- summary.grid[i,"meanBkm2"]

    if(all(is.na(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]])) ){
      cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- "#FFFFFF"
    } else {

      if(all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[2] ) & all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] >= cat_lim[1])){
        cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[1]
      }

      if(all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[3])  & all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[2])){
        cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[2]
      }

      if(all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[4])  & all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[3])){
        cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[3]
      }

      if(all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[5])  & all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[4])){
        cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[4]
      }

      if(all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[6])  & all(cgpmgrid$meanBkm2[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[5])){
        cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[5]
      }
    }
  }


  lx =map_range["xmax"] - map_range["xmin"]
  ly =map_range["ymax"] - map_range["ymin"]
  ratio <- ly/lx
  img_width <- 30
  img_height <- img_width*ratio

  # par( mar=c(4, 5, 4, 2)) #c(bottom, left, top, right)
  # plot(1,1,type="n",xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), xlab=expression(paste("Longitude (",degree,"E)")), ylab=expression(paste("Latitude (",degree,"N)")), main=paste("Biomass (kg/km^2) -", sspp), asp=1 )
  # plot(cgpmgrid,col=cgpmgrid$color, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")                                                                                 #data = centroidi_coords,
  # if (stratum == "200,800" | stratum == 3){
  #   plot(mask,col="white", xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")
  # }
  # plot(continent, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), border="grey", col="light grey", add=T)
  #
  # legend("topleft", cate_label, col=palette_cate, pch = 15, pt.cex=c(2), title = "Biomass (kg/km^2)\n(classes based on percentiles)", bty = "n")
  # box()
  # p3 <- recordPlot()

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
  # cgpmgrid$level <- as.character( cgpmgrid$level)
  cgpmgrid$level <- is.factor(cgpmgrid$level)
  cgpmgrid$level <- factor(cgpmgrid$level, levels=seq(1,nrow(tab_color),1))
  cgpmgrid <- cgpmgrid %>% arrange(level)

  p3 <- ggplot(cgpmgrid)+
    geom_spatvector(aes(fill = category))+
    scale_fill_manual(values = cols, breaks=names(cols))+
    geom_polygon(data=world, aes(long,lat,group=group), fill="dark grey", colour="darkgrey")+
    coord_sf(xlim = xl, ylim = yl)+
    labs(fill="Biomass\n (kg/km^2)")+
    xlab("longitude")+
    ylab("latitude")

  print(p3)

  if(save){
    jpeg(filename=paste(wd, "/output/",sspp," - GFCM GRID BIOMASS.jpg", sep = ""), width=25, height=25, bg="white", units="cm",res=200)
    print(p3)
    dev.off()
  }



  ###################
  ###################
  ###             ###
  ###     MIW     ###
  ###             ###
  ###################
  ###################

  cgpmgrid <- cgpmgrid_bkp
  summary.grid <- meanWEIGHT_grid
  # metaDB <- mTATB

  #----------------------------------------
  # Saving maps of the MIW
  #----------------------------------------                                                                          # [summary.grid$meanNkm2>5]
  cat_lim <- as.numeric(quantile(summary.grid[!is.na(summary.grid$MIW), "MIW"], probs = seq(0, 1,0.20) ) )
  #cat_lim[1] <- 0
  #  cate_label[len] <- ifelse( round(cat_lim[len],0) != 0, round(cat_lim[len],0), ifelse( round(cat_lim[len],1) != 0,  round(cat_lim[len],1),  round(cat_lim[len],2)))
  cate_label <-   round(cat_lim,3)
  cate_label_min <- cate_label
  #cate_label_min[1] <- 0
  cate_label_max <- cate_label[2:length(cate_label)]


  cate_label <- paste("(", cate_label_min, ", ", cate_label_max, "]", sep="")
  cate_label[1] <- paste("[", cate_label_min[1], ", ", cate_label_max[1], "]", sep="")

  #cate_label[length(cate_label)] <- paste(">", cate_label_max[length(cate_label_max)] )
  cate_label <- cate_label[1: (length(cate_label)-1)]

  palette_cate <- heat.colors(length(cate_label))
  palette_cate <- palette_cate[length(palette_cate):1]

  cgpmgrid$MIW <- NA
  cgpmgrid$color <- "#FFFFFF"

  i=3
  for (i in 1:length(summary.grid[,1])){
    cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- summary.grid[i,"MIW"]

    if(all(is.na(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]])) ){
      cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- "#FFFFFF"
    } else {

      if(all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[2])  & all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] >= cat_lim[1])){
        cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[1]
      }

      if(all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[3])  & all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[2])){
        cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[2]
      }

      if(all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[4])  & all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[3])){
        cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[3]
      }

      if(all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[5])  & all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[4])){
        cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[4]
      }

      if(all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <= cat_lim[6])  & all(cgpmgrid$MIW[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] > cat_lim[5])){
        cgpmgrid$color[cgpmgrid$GFCM_ID==summary.grid[i,"cgpmgridlevel"]] <- palette_cate[5]
      }
    }
  }


  lx =map_range["xmax"] - map_range["xmin"]
  ly =map_range["ymax"] - map_range["ymin"]
  ratio <- ly/lx
  img_width <- 30
  img_height <- img_width*ratio


  # par( mar=c(4, 5, 4, 2)) #c(bottom, left, top, right)
  # plot(1,1,type="n",xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]) ,  xlab=expression(paste("Longitude (",degree,"E)")), ylab=expression(paste("Latitude (",degree,"N)")), main=paste(paste("Mean Individual Weight  (kg)", sep=" "), sspp), asp=1)
  # plot(cgpmgrid, col=cgpmgrid$color, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")
  # if (stratum == "200,800" | stratum == 3){
  #   plot(mask,col="white", xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), add=T, border="light grey")
  # }
  # plot(continent, xlim=c(map_range["xmin"],map_range["xmax"]), ylim=c(map_range["ymin"],map_range["ymax"]), border="grey", col="light grey", add=T)
  #
  # legend("topleft", cate_label, col=palette_cate, pch = 15, pt.cex=c(2), title = "MIW (kg)\n(classes based on percentiles)", bty = "n")
  #
  # box()
  # p4 <- recordPlot()

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
  # cgpmgrid$level <- as.character( cgpmgrid$level)
  cgpmgrid$level <- is.factor(cgpmgrid$level)
  cgpmgrid$level <- factor(cgpmgrid$level, levels=seq(1,nrow(tab_color),1))
  cgpmgrid <- cgpmgrid %>% arrange(level)

  p4 <- ggplot(cgpmgrid)+
    geom_spatvector(aes(fill = category))+
    scale_fill_manual(values = cols, breaks=names(cols))+
    geom_polygon(data=world, aes(long,lat,group=group), fill="dark grey", colour="darkgrey")+
    coord_sf(xlim = xl, ylim = yl)+
    labs(fill="MIW (kg)")+
    xlab("longitude")+
    ylab("latitude")

  print(p4)

  if(save){
    jpeg(filename=paste(wd, "/output/",sspp," - GFCM GRID MIW.jpg", sep = ""), width=25, height=25, bg="white", units="cm",res=200)
    print(p4)
    dev.off()
  }


  return(list(abundance_grid, biomass_grid, meanWEIGHT_grid))
}

