#' Bubble plot of abundance and biomass indices by haul
#'
#' @author Walter Zupa \email{zupa@fondazionecoispa.org}
#' @description
#' The function generates bubble plot of abundance and biomass indices by haul.
#' If no resolution is specified (res = NA), the function works offline and uses an internal bathymetry dataset (med_bathy) covering the Mediterranean and Black Sea, reducing the computational time.
#'
#' @param mTATB data frame
#' @param map_lim coordinates limits for the plotted map
#' @param depth_lines vector of three depth bathymetrical lines to be plotted
#' @param buffer buffer to the coordinate limits in map units
#' @param res resolution of the bathymetrical lines
#' @param wd working directory
#' @param save boolean. If TRUE the plot is saved in the user defined working directory (wd)
#' @param verbose boolean. If TRUE a message is printed
#' @importFrom marmap getNOAA.bathy as.xyz
#' @importFrom utils data
#' @importFrom ggplot2 coord_sf geom_polygon scale_x_continuous scale_y_continuous geom_contour geom_point scale_size coord_map ggtitle theme ggsave element_blank element_rect element_text map_data aes labs
#' @return A \code{ggplot} object representing the bubble plot of abundance or biomass by hauls.
#' @examples
#' \donttest{
#' # Use internal data
#' data(TA)
#' data(TB)
#' # Create a merged dataset for GSA 10 and year 2016
#' mTATB <- merge_TATB(TA[TA$AREA == 10 & TA$YEAR == 2016, ],
#'                     TB[TB$AREA == 10 & TB$YEAR == 2016, ],
#'                     "MERLMER")
#' bubble_plot_by_haul_indexes(mTATB, map_lim = NA, depth_lines = c(200, 500, 800), save = FALSE)
#' }
#' @export

bubble_plot_by_haul_indexes <- function(mTATB, map_lim, depth_lines, buffer = 0, res = NA, wd = NA, save = TRUE, verbose = TRUE) {
  if (is.na(wd)) {
    wd <- tempdir()
  }

  if (FALSE) {
    map_lim <- c(15.5, 20.0, 39.8, 42.5)
    # map_lim <- NA
    depth_lines <- c(200, 500, 800)
    buffer <- 0.5
    res <- NA
    save <- TRUE
    verbose <- TRUE
    # ta <- read.table(paste(wd,"input/TA/TA GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
    # tb <- read.table(paste(wd,"input/TB/TB GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
    # tc <- read.table(paste(wd,"input/TC/TC GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
    wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\R_BioIndex_3.3_(in update)"

    mTATB <- read.table("D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\R_BioIndex_3.3_(in update)\\output\\mergeTATB_MERLMER.csv", sep = ";", header = TRUE)

    bubble_plot_by_haul_indexes(mTATB, map_lim, depth_lines, buffer = 0, wd = wd, save = TRUE, verbose = TRUE)
  }

  if (!requireNamespace("mapproj", quietly = TRUE)) {
    stop("Package 'mapproj' needed for coord_map() in ggplot2. Please install it.", call. = FALSE)
  }

  if (is.na(wd) & save) {
    if (verbose) {
      message("Missing working directory. Results are not saved in the local folder.")
    }
  }

  N_km2 <- V1 <- V2 <- V3 <- group <- kg_km2 <- lat <- lon <- long <- NULL

  dep_text_N <- expression(paste("Abundance by haul ", (n / km^2), sep = " ")) # ," - ",depth_range[1],"-",depth_range[2]," m"
  dep_text_kg <- expression(paste("Biomass by haul ", (kg / km^2), sep = " ")) # ," - ",depth_range[1],"-",depth_range[2]," m"

  metaDB <- mTATB
  GENERE <- as.character(unique(metaDB$GENUS)[!is.na(unique(metaDB$GENUS) != -1) & unique(metaDB$GENUS) != -1])
  SPECIE <- as.character(unique(metaDB$SPECIES)[!is.na(unique(metaDB$SPECIES) != -1) & unique(metaDB$SPECIES) != -1])
  sspp <- paste(GENERE, SPECIE, sep = "")
  GSA <- unique(metaDB$GSA)
  species <- sspp


  if (any(is.na(map_lim))) {
    map_lim[c(1:2)] <- range(metaDB$MEAN_LONGITUDE_DEC, na.rm = TRUE)
    map_lim[c(3, 4)] <- range(metaDB$MEAN_LATITUDE_DEC, na.rm = TRUE)
    lx <- map_lim[2] - map_lim[1]
    ly <- map_lim[4] - map_lim[3]
  } else {
    lx <- map_lim[2] - map_lim[1]
    ly <- map_lim[4] - map_lim[3]
  }

  x1 <- map_lim[1]
  x2 <- map_lim[2]
  y1 <- map_lim[3]
  y2 <- map_lim[4]

  map_ratio <- ly / lx
  theme_opts <- list(theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "light blue", linetype = "solid", color = "black"),
    plot.background = element_rect(
      fill = "white",
      linewidth = 1, linetype = "solid", color = "black"
    ),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 22)
  ))


  buff <- buffer


  world <- map_data("world")
  xl <- c(x1, x2)
  yl <- c(y1, y2)
  x_breaks <- c(round(xl[1], 0), round(xl[1], 0) + round((xl[2] - xl[1]) / 2, 0), round(xl[1], 0) + 2 * round((xl[2] - xl[1]) / 2, 0))
  y_breaks <- c(round(yl[1], 0), round(yl[1], 0) + round((yl[2] - yl[1]) / 2, 0), round(yl[1], 0) + 2 * round((yl[2] - yl[1]) / 2, 0))


  # bath <- suppressMessages(getNOAA.bathy (lon1 = min(xl)-1, lon2 = max(xl)+1, lat1 = min(yl)-1-buff,lat2 = max(yl)+1+buff, resolution = res))
  # bat_xyz <- as.xyz(bath)

  if (!is.na(res)) {
    bath <- suppressMessages(getNOAA.bathy(lon1 = min(xl) - 1, lon2 = max(xl) + 1, lat1 = min(yl) - 1 - buff, lat2 = max(yl) + 1 + buff, resolution = res))
    bat_xyz <- as.xyz(bath)
  } else {
    data("med_bathy", package = "BioIndex", envir = environment())
    bat_xyz <- as.xyz(med_bathy)
    bat_xyz <- bat_xyz[!is.na(bat_xyz[, 3]), ]
    colnames(bat_xyz) <- c("V1", "V2", "V3")
    bat_xyz <- subset(bat_xyz, V1 >= min(xl) & V1 <= max(xl) & V2 >= min(yl) - buff & V2 <= max(yl) + buff)
  }

  if (any(is.na(depth_lines))) {
    depth_lines <- c(-200, -500, -800)
  } else {
    depth_lines <- -1 * depth_lines
  }


  if (nrow(metaDB) != 0) {
    loc <- data.frame(id = metaDB$id, COUNTRY = metaDB$COUNTRY, GSA = metaDB$GSA, YEAR = metaDB$YEAR, lon = metaDB$MEAN_LONGITUDE_DEC, lat = metaDB$MEAN_LATITUDE_DEC, N_km2 = metaDB$N_km2, kg_km2 = metaDB$kg_km2)

    loc1 <- loc[loc$N_km2 != 0, ]
    br_n <- c(0, as.numeric(summary(loc1$N_km2)[-4]))
    br_kg <- c(0, as.numeric(summary(loc1$kg_km2)[-4]))
    labels_n <- c("0", as.character(round(as.numeric(summary(loc1$N_km2))[-4], 3)))
    labels_kg <- c("0", as.character(round(as.numeric(summary(loc1$kg_km2))[-4], 3)))

    # coordinates(loc)<-c("lon","lat")
    # proj4string(loc) <- CRS("+proj=longlat")
    # loc_laea<-spTransform(loc, CRS("+proj=longlat"))
    # loc_laea_df<-data.frame(loc_laea)
    loc[loc$N_km2 == 0, "colore"] <- "red"
    loc[loc$N_km2 != 0, "colore"] <- "blue"

    limits <- data.frame(matrix(NA, ncol = 2, nrow = 2))
    colnames(limits) <- c("X", "Y")
    limits[1, 1] <- x1
    limits[1, 2] <- y1
    limits[2, 1] <- x2
    limits[2, 2] <- y2
    # coordinates(limits) <- ~ X+Y
    # proj4string(limits) <- CRS("+proj=longlat")
    # lim_laea<-as.data.frame(spTransform(limits, CRS("+proj=longlat")))
    xmin <- x1
    xmax <- x2
    ymin <- y1
    ymax <- y2


    # ABUNDANCE OUTPUT
    suppressWarnings(suppressMessages(
      p_n <- ggplot() +
        # coord_map() requires the mapproj package to be installed.
        coord_map(projection = "mercator", xlim = c((xmin - buff), (xmax + buff)), ylim = c((ymin - buff), (ymax + buff))) +
        geom_polygon(data = world, aes(long, lat, group = group), fill = "light grey", colour = "darkgrey") +
        geom_contour(
          data = bat_xyz,
          aes(x = V1, y = V2, z = V3),
          breaks = depth_lines[3], color = "#1F618D", linewidth = 0.5
        ) +
        geom_contour(
          data = bat_xyz,
          aes(x = V1, y = V2, z = V3),
          breaks = depth_lines[2], color = "#5499C7", linewidth = 0.5
        ) +
        geom_contour(
          data = bat_xyz,
          aes(x = V1, y = V2, z = V3),
          breaks = depth_lines[1], color = "#4db5fa", linewidth = 0.5
        ) +
        geom_point(
          data = loc, aes(lon, lat, group = NA, fill = NA, size = N_km2), shape = 21,
          colour = "black", fill = loc$colore, alpha = I(3 / 10)
        ) +
        scale_size(range = c(1, 14), breaks = br_n, labels = labels_n, guide = "legend", labs(size = "n/km^2")) +
        # theme(aspect.ratio=map_ratio*1.5)+
        ggtitle(dep_text_N) +
        theme_opts
    ))

    # print(p_n) # rimosso per policy CRAN
    if (verbose) {
      message("Plot of abundance by haul generated successfully.")
    }
    if (save) {
      if (is.na(wd)) {
        warning("\nNo working directory defined by the user. The plot was not saved in the local folder\n")
      } else if (!is.na(wd)) {
        ggsave(paste(wd, "/output/", sspp, "_GSA", GSA, "-Abundance_by_haul.jpg", sep = ""),
          plot = p_n, width = 20, height = 20, units = "cm", dpi = 300
        )
        if (verbose) {
          message("Plot of abundance by haul saved successfully")
        }
      }
    }

    # BIOMASS OUTPUT
    suppressWarnings(suppressMessages(
      p_k <- ggplot() +
        coord_map(projection = "mercator", xlim = c((xmin - buff), (xmax + buff)), ylim = c((ymin - buff), (ymax + buff))) +
        geom_polygon(data = world, aes(long, lat, group = group), fill = "light grey", colour = "darkgrey") +
        geom_contour(
          data = bat_xyz,
          aes(x = V1, y = V2, z = V3),
          breaks = depth_lines[3], color = "#1F618D", linewidth = 0.5
        ) +
        geom_contour(
          data = bat_xyz,
          aes(x = V1, y = V2, z = V3),
          breaks = depth_lines[2], color = "#5499C7", linewidth = 0.5
        ) +
        geom_contour(
          data = bat_xyz,
          aes(x = V1, y = V2, z = V3),
          breaks = depth_lines[1], color = "#4db5fa", linewidth = 0.5
        ) +
        geom_point(
          data = loc, aes(lon, lat, group = NA, fill = NA, size = kg_km2), shape = 21,
          colour = "black", fill = loc$colore, alpha = I(3 / 10)
        ) +
        scale_size(range = c(1, 14), breaks = br_n, labels = labels_n, guide = "legend", labs(size = "kg/km^2")) +
        # theme(aspect.ratio=map_ratio*1.5)+
        ggtitle(dep_text_kg) +
        theme_opts
    ))

    # print(p_k) # rimosso per policy CRAN
    if (verbose) {
      message("Plot of biomass by haul generated successfully.")
    }
    if (save) {
      if (is.na(wd)) {
        warning("\nNo working directory defined by the user. The plot was not saved in the local folder\n")
      } else if (!is.na(wd)) {
        ggsave(paste(wd, "/output/", sspp, "_GSA", GSA, "-Biomass_by_haul.jpg", sep = ""),
          plot = p_k, width = 20, height = 20, units = "cm", dpi = 300
        )
        if (verbose) {
          message("Plot of biomass by haul saved successfully")
        }
      }
    }

    return(list(p_n, p_k))
  } else {
    if (verbose) {
      message("Not enough data to be plotted\n")
      message("Bubble plots - indices  by haul skipped\n\n")
    }
  }
}
