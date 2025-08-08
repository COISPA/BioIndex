#' Plot of hauls time series
#'
#' @description
#' This function generates a spatial plot of haul positions over time, displaying haul labels and bathymetric lines over a customizable map extent.
#' If no resolution is specified (res = NA), the function works offline and uses an internal bathymetry dataset (med_bathy) covering the Mediterranean and Black Sea, reducing the computational time.
#' @param mTATB data frame
#' @param country country code as reported in MEDITS format. "all" code to perform the analysis on all the countries of the same GSA
#' @param map_lim coordinates limits for the plotted map
#' @param depth_lines vactor of three depth bathymetrical lines to be plotted
#' @param buffer buffer to the coordinate limits in map units
#' @param res resolution of the bathymetrical lines
#' @param wd working directory
#' @param save boolean. If TRUE the plot is saved in the user defined working directory (wd)
#' @param verbose boolean. If TRUE a message is printed
#' @importFrom marmap getNOAA.bathy as.xyz
#' @importFrom ggplot2 ggplot coord_sf geom_polygon scale_x_continuous scale_y_continuous geom_contour geom_text geom_point coord_map ggtitle ggsave theme map_data
#' @export
hauls_position <- function(mTATB, country="all",map_lim,depth_lines, buffer=0, res=NA, wd=NA,save=TRUE, verbose=TRUE){

  Haul <- V1 <- V2 <- V3 <- Year <- group <- merge_TATB <-  lat <- lon <- long <-NULL #

  if (FALSE) {
    map_lim <- c(15.5,20.0,39.8,42.5)
    map_lim <- NA
    depth_lines <- c(200,500,800)
    buffer=0.5
    res=0.1
    save=TRUE
    ta <- read.table(paste(wd,"input/TA/TA GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
    tb <- read.table(paste(wd,"input/TB/TB GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
    tc <- read.table(paste(wd,"input/TC/TC GSA18 2017-2020.csv", sep="/"), sep=";", header=T)
    wd <- "D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\R_BioIndex_3.3_(in update)"

    m <- read.table("D:\\Documents and Settings\\Utente\\Documenti\\GitHub\\BioIndex\\R_BioIndex_3.3_(in update)\\output\\mergeTATB_MERLMER.csv",sep=";",header=TRUE)

    hauls_position(mTATB=m,map_lim,depth_lines, buffer=0.1, wd=wd,save=FALSE, verbose=TRUE)
  }

  if (!requireNamespace("mapproj", quietly = TRUE)) {
    stop("Package 'mapproj' needed for this function to work. Please install it.", call. = FALSE)
  }

  if (is.na(wd) & save) {
    save =FALSE
    if (verbose){
      message("Missing working directory. Results are not saved in the local folder.")
    }
  }




countries <- unique(mTATB[!is.na(mTATB$COUNTRY),"COUNTRY"])
l_country <- length(unique(mTATB$COUNTRY))

if (l_country==1){
  check_country ="Y"
  country_analysis <- countries
  if (verbose){
    message(paste0("Analysis conducted for the following country: ",countries))
  }
} else {
  if (any(country %in% "all")) {
    country_analysis <<- countries
    if (verbose){
      message(paste0("Analysis conducted for the following countries: ",paste(countries,collapse = ", ")))
    }
  } else {
    country_analysis <<- country
    if (verbose){
      message(paste0("Analysis conducted for the following countries: ",paste(countries,collapse = ", ")))
    }
  }
}

mTATB <- mTATB[mTATB$COUNTRY %in% country_analysis,  ]

metaDB <- mTATB

GENERE <- as.character(unique(metaDB$GENUS)[ !is.na(unique(metaDB$GENUS) != -1) & unique(metaDB$GENUS) != -1])
SPECIE <- as.character(unique(metaDB$SPECIES)[!is.na(unique(metaDB$SPECIES) != -1) & unique(metaDB$SPECIES) != -1])
sspp <- paste(GENERE,SPECIE, sep="")
GSA <- unique(metaDB$GSA)
species <- sspp
#### selection of map limits by file ####
      if (any(is.na(map_lim))) {
        map_lim[c(1:2)] <- range(metaDB$MEAN_LONGITUDE_DEC,na.rm=TRUE)
        map_lim[c(3,4)] <- range(metaDB$MEAN_LATITUDE_DEC,na.rm=TRUE)
        lx =map_lim[2] - map_lim[1]
        ly =map_lim[4] - map_lim[3]
      } else {
        lx =map_lim[2] - map_lim[1]
        ly =map_lim[4] - map_lim[3]
      }


x1 <- map_lim[1]
x2 <- map_lim[2]
y1 <- map_lim[3]
y2 <- map_lim[4]

map_ratio <- ly/lx

ddd <- metaDB

loc <- data.frame(id=ddd$id, Haul = ddd$HAUL_NUMBER, Year = as.factor(ddd$YEAR), lon=ddd$MEAN_LONGITUDE_DEC, lat=ddd$MEAN_LATITUDE_DEC)

theme_opts<-list(theme(panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.background = element_rect(fill = 'light blue',linetype="solid",color="black"),
                       plot.background = element_rect(fill="white",
                                                      linewidth=1,linetype="solid",color="black"),
                       axis.line = element_blank(),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks = element_blank(),
                       axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       plot.title = element_text(size=22)))


buff<-buffer
dep_text <-"Hauls position"
world <- map_data("world")
xl <- c(x1,x2)
yl <- c(y1,y2)
x_breaks <- c(round(xl[1], 0), round(xl[1], 0) + round((xl[2] - xl[1]) / 2, 0), round(xl[1], 0) + 2 * round((xl[2] - xl[1]) / 2, 0))
y_breaks <- c(round(yl[1], 0), round(yl[1], 0) + round((yl[2] - yl[1]) / 2, 0), round(yl[1], 0) + 2 * round((yl[2] - yl[1]) / 2, 0))

lon_1 = min(xl)-1-buff
lon_2 = max(xl)+1+buff
lat_1 = min(yl)-1-buff
lat_2 = max(yl)+1+buff

if (!is.na(res)){
  bath <- suppressMessages(getNOAA.bathy(lon1 = min(xl) - 1, lon2 = max(xl) + 1, lat1 = min(yl) - 1 - buff, lat2 = max(yl) + 1 + buff, resolution = res))
  bat_xyz <- as.xyz(bath)
} else {
  data("med_bathy", package = "BioIndex", envir = environment())
  bat_xyz <- as.xyz(med_bathy)
  bat_xyz <- bat_xyz[!is.na(bat_xyz[,3]), ]
  colnames(bat_xyz) <- c("V1","V2","V3")
  bat_xyz <- subset(bat_xyz, V1 >= min(xl) & V1 <= max(xl) & V2 >= min(yl) - buff & V2 <= max(yl) + buff)
}


if (any(is.na(depth_lines))) {
  depth_lines <- c(-200,-500,-800)
} else {
  depth_lines <- -1 * depth_lines
}


if (any(is.na(depth_lines))) {
  depth_lines <- c(-200,-500,-800)
} else {
  depth_lines <- -1 * depth_lines
}

p <- suppressMessages(
  ggplot() +
  coord_sf(xlim = xl, ylim = yl, expand = TRUE) +
  geom_polygon(data=world, aes(long,lat,group=group), fill="grey")+
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks) +
  geom_contour(data = bat_xyz,
               aes(x = V1, y = V2, z = V3),
               breaks = depth_lines[3], color = "#1F618D", linewidth = 0.5) +
  geom_contour(data = bat_xyz,
               aes(x = V1, y = V2, z = V3),
               breaks = depth_lines[2], color = "#5499C7", linewidth = 0.5) +
  geom_contour(data = bat_xyz,
               aes(x = V1, y = V2, z = V3),
               breaks = depth_lines[1], color = "#4db5fa", linewidth = 0.5) +
  geom_text(data=loc, aes(lon, lat, label=Haul)) +
  geom_point(data=loc, aes(lon, lat,color=Year, fill=Year)) +
  coord_map(projection="mercator",xlim = c((x1-buff),(x2+buff)), ylim = c((y1-buff),(y2+buff))) +
  ggtitle(dep_text) +
  theme_opts)
print(p)
if (verbose) {
  message("Plot of Hauls position generated successfully.")
}


if (save){
  if (is.na(wd)){
    warning("\nNo working directory defined by the user. The plot was not saved in the local folder\n")
  } else if (!is.na(wd)) {
    ggsave(paste(wd, "/output/",sspp,"_GSA",GSA," -Hauls position.jpg", sep=""),
           width = 20, height = 20, units ="cm", dpi = 300)
    if(verbose){
      message("Bubble plot of Hauls position correctly saved \n")
    }
  }
}

return(p)
}
