#' Bubbleplot of abundance indices for recruits and spawners
#'
#' @author Walter Zupa \email{zupa@fondazionecoispa.org}
#' @description
#' The function generates bubbleplots of abundance indices for recruits and spawners
#' If no resolution is specified (res = NA), the function works offline and uses an internal bathymetry dataset (med_bathy) covering the Mediterranean and Black Sea, reducing the computational time.
#'
#' @param mTATC mTATC table
#' @param map_range range of coordinates for the map
#' @param thresh_rec threshold value to select recruits data from mTATC table (reported in mm)
#' @param thresh_spaw threshold value to select spawners data from mTATC table (reported in mm)
#' @param depths three reference bathymetric lines to be plotted in the maps
#' @param res resolution of the depth lines
#' @param buffer buffer around the map
#' @param wd working directory
#' @param save boolean. If TRUE the outputs are saved in the local folder
#' @param verbose boolean. If TRUE messages are prompted in the console
#' @return A \code{ggplot} object representing the bubble plot of sex ratio by hauls.
#' @examples
#' \donttest{
#' # Create a merged dataset for GSA 10
#' mTATC <- merge_TATC(TA[TA$AREA==10,], TC[TC$AREA==10,], "MERLMER")
#' bubbleplot_RS_by_hauls(mTATC, map_range=c(15,21,39,43), thresh_rec=200, thresh_spaw=250)
#' }
#' @export
#' @importFrom marmap getNOAA.bathy as.xyz
#' @importFrom utils data
#' @importFrom ggplot2 coord_sf geom_polygon scale_x_continuous scale_y_continuous geom_contour geom_point scale_size coord_map ggtitle theme ggsave element_blank element_rect element_text map_data aes labs
#' @importFrom stats aggregate
bubbleplot_RS_by_hauls <- function(mTATC, map_range, thresh_rec, thresh_spaw, depths = c(50, 200, 800), res=NA, buffer=0.1,wd=NA, save=TRUE, verbose = FALSE) {
    if (is.na(wd)) { wd <- tempdir() }
  if (FALSE) {
    thresh_rec <- 20
    thresh_spaw <- 21
    map_range <- c(15.0, 21.0, 39.5, 42.5)
    depths = c(50, 200, 800)
    res=1
    buffer=0.1
    save <- TRUE
    verbose <- TRUE
    bubbleplot_RS_by_hauls(mTATC, map_range, thresh_rec, thresh_spaw, depths = c(50, 200, 800))
  }

  if (!requireNamespace("mapproj", quietly = TRUE)) {
    stop("Package 'mapproj' needed for this function to work. Please install it.", call. = FALSE)
  }

  long <- lat <- group <- V1 <- V2 <- V3 <- lon <- indices <- NULL

  buff <- buffer

  merge_TATC <- mTATC
  GENERE <- as.character(unique(merge_TATC$GENUS)[unique(merge_TATC$GENUS) != -1])
  SPECIE <- as.character(unique(merge_TATC$SPECIES)[unique(merge_TATC$SPECIES) != -1])
  sspp <- paste(GENERE, SPECIE, sep = "")

  depths <- -1 * depths


  names(map_range) <- c("xmin", "xmax", "ymin", "ymax")
  x1 <- min(map_range[1])
  x2 <- max(map_range[2])
  y1 <- min(map_range[3])
  y2 <- max(map_range[4])

  xl <- c(x1, x2)
  yl <- c(y1, y2)

  map_lim <- map_range
  lx <- map_range["xmax"] - map_range["xmin"]
  ly <- map_range["ymax"] - map_range["ymin"]
  map_ratio <- ly / lx

  limits <- data.frame(matrix(NA, ncol = 2, nrow = 2))
  colnames(limits) <- c("X", "Y")
  limits[1, 1] <- map_lim[1]
  limits[1, 2] <- map_lim[3]
  limits[2, 1] <- map_lim[2]
  limits[2, 2] <- map_lim[4]

  GSA <- unique(merge_TATC$GSA)[[1]]
  threshold <- thresh_rec

  ddd_r <- NA
  ddd_r <- merge_TATC[merge_TATC$LENGTH_CLASS <= threshold & merge_TATC$LENGTH_CLASS != -1, ]

  # load("data/med_bathy.rda")
  if (!is.na(res)){
    bath <- suppressMessages(suppressWarnings(getNOAA.bathy(lon1 = min(xl) - 1, lon2 = max(xl) + 1, lat1 = min(yl) - 1 - buff, lat2 = max(yl) + 1 + buff, resolution = res)))
    bat_xyz <- as.xyz(bath)
  } else {
  data("med_bathy", package = "BioIndex", envir = environment())
  bat_xyz <- as.xyz(med_bathy)
  bat_xyz <- bat_xyz[!is.na(bat_xyz[,3]), ]
  colnames(bat_xyz) <- c("V1","V2","V3")
  bat_xyz <- subset(bat_xyz, V1 >= min(xl) & V1 <= max(xl) & V2 >= min(yl) - buff & V2 <= max(yl) + buff)
  }

  if (nrow(ddd_r) != 0) {
    pivot_r <- aggregate(as.numeric(as.character(ddd_r$N_km2)), by = list(ddd_r$id, ddd_r$MEAN_LONGITUDE_DEC, ddd_r$MEAN_LATITUDE_DEC), sum)
    colnames(pivot_r) <- c("id", "lon", "lat", "indices")
    loc <- pivot_r
    br <- as.numeric(summary(loc$indices)[-4])
    labels <- as.character(round(as.numeric(summary(loc$indices))[-4], 3))

    theme_opts <- list(theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "light blue", linetype = "solid", color = "black"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 22)
    ))


    dep_text <- expression(paste("Abundance of recruits ", (n / km^2), sep = " "))

    world <- map_data("world")
    x_breaks <- c(round(xl[1], 0), round(xl[1], 0) + round((xl[2] - xl[1]) / 2, 0), round(xl[1], 0) + 2 * round((xl[2] - xl[1]) / 2, 0))
    y_breaks <- c(round(yl[1], 0), round(yl[1], 0) + round((yl[2] - yl[1]) / 2, 0), round(yl[1], 0) + 2 * round((yl[2] - yl[1]) / 2, 0))


    suppressWarnings(suppressMessages(
      pr <- ggplot() +
        coord_sf(xlim = xl, ylim = yl, expand = TRUE) +
        geom_polygon(data = world, aes(long, lat, group = group), fill = "grey") +
        scale_x_continuous(breaks = x_breaks) +
        scale_y_continuous(breaks = y_breaks) +
        geom_contour(
          data = bat_xyz,
          aes(x = V1, y = V2, z = V3),
          breaks = depths[3], color = "#1F618D", linewidth = 0.5
        ) +
        geom_contour(
          data = bat_xyz,
          aes(x = V1, y = V2, z = V3),
          breaks = depths[2], color = "#5499C7", linewidth = 0.5
        ) +
        geom_contour(
          data = bat_xyz,
          aes(x = V1, y = V2, z = V3),
          breaks = depths[1], color = "#4db5fa", linewidth = 0.5
        ) +
        geom_point(
          data = loc, aes(lon, lat, group = NULL, fill = NULL, size = indices), shape = 21,
          color = "black", fill = "blue", alpha = I(4 / 10)
        ) +
        scale_size(range = c(1, 14), breaks = br, labels = labels, guide = "legend", labs(size = "n/km^2")) +
        coord_map(projection = "mercator", xlim = c((x1 - buff), (x2 + buff)), ylim = c((y1 - buff), (y2 + buff))) +
        ggtitle(dep_text) +
        theme_opts
    ))
        # print(pr) # rimosso per policy CRAN
    if (save) {
      ggsave(paste(wd, "/output/", sspp, "_GSA", GSA, " -indices of RECRUITS.jpg", sep = ""),
        width = 20, height = 20, units = "cm", dpi = 300, plot = pr
      )
    }
    if (verbose) {
      message("Bubble plot of recruits correctly saved ")
    }
  } else {
    if (verbose) {
      message("Not enough data for recruits to be plotted")
      message("Bubble plots - indices of recruits skipped\n")
    }
  }

  ##################
  ##################
  ##################

  merge_TATC <- mTATC
  threshold <- thresh_spaw
  ddd_s <- NA
  ddd_s <- merge_TATC[merge_TATC$LENGTH_CLASS >= threshold & merge_TATC$LENGTH_CLASS != -1 & merge_TATC$SEX == "F", ]

  if (nrow(ddd_s) != 0) {
    pivot_s <- aggregate(as.numeric(as.character(ddd_s$N_km2)), by = list(ddd_s$id, ddd_s$MEAN_LONGITUDE_DEC, ddd_s$MEAN_LATITUDE_DEC), sum)
    colnames(pivot_s) <- c("id", "lon", "lat", "indices")
    loc <- pivot_s
    br <- as.numeric(summary(loc$indices)[-4])
    labels <- as.character(round(as.numeric(summary(loc$indices)[-4]), 3))

    theme_opts <- list(theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "light blue", linetype = "solid", color = "black"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 22)
    ))

    dep_text <- expression(paste("Abundance of spawners ", (n / km^2), sep = " "))

    # res <- 0.1
    world <- map_data("world")
    x_breaks <- c(round(xl[1], 0), round(xl[1], 0) + round((xl[2] - xl[1]) / 2, 0), round(xl[1], 0) + 2 * round((xl[2] - xl[1]) / 2, 0))
    y_breaks <- c(round(yl[1], 0), round(yl[1], 0) + round((yl[2] - yl[1]) / 2, 0), round(yl[1], 0) + 2 * round((yl[2] - yl[1]) / 2, 0))

    # bath <- suppressMessages(getNOAA.bathy(lon1 = min(xl) - 1, lon2 = max(xl) + 1, lat1 = min(yl) - 1 - buff, lat2 = max(yl) + 1 + buff, resolution = res))
    # bat_xyz <- as.xyz(bath)

    suppressWarnings(suppressMessages(
      ps <- ggplot() +
        coord_sf(xlim = xl, ylim = yl, expand = TRUE) +
        geom_polygon(data = world, aes(long, lat, group = group), fill = "grey") +
        scale_x_continuous(breaks = x_breaks) +
        scale_y_continuous(breaks = y_breaks) +
        geom_contour(
          data = bat_xyz,
          aes(x = V1, y = V2, z = V3),
          breaks = depths[3], color = "#1F618D", linewidth = 0.5
        ) +
        geom_contour(
          data = bat_xyz,
          aes(x = V1, y = V2, z = V3),
          breaks = depths[2], color = "#5499C7", linewidth = 0.5
        ) +
        geom_contour(
          data = bat_xyz,
          aes(x = V1, y = V2, z = V3),
          breaks = depths[1], color = "#4db5fa", linewidth = 0.5
        ) +
        geom_point(
          data = loc, aes(lon, lat, group = NULL, fill = NULL, size = indices), shape = 21,
          color = "black", fill = "blue", alpha = I(4 / 10)
        ) +
        scale_size(range = c(1, 14), breaks = br, labels = labels, guide = "legend", labs(size = "n/km^2")) +
        coord_map(projection = "mercator", xlim = c((x1 - buff), (x2 + buff)), ylim = c((y1 - buff), (y2 + buff))) +
        ggtitle(dep_text) +
        theme_opts
    ))
        # print(ps) # rimosso per policy CRAN
    if (save) {
      ggsave(paste(wd, "/output/", sspp, "_GSA", GSA, " -indices of SPAWNERS.jpg", sep = ""),
        width = 20, height = 20, units = "cm", dpi = 300, plot = ps
      )
    }
    if (verbose) {
      message("Bubble plot of spawners correctly saved ")
    }
  } else {
    if (verbose) {
      message("Not enough data for spawners to be plotted")
      message("Bubble plots - indices of spawners skipped")
    }
  }
}