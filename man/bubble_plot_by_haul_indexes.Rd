% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bubble_plot_by_haul_indexes.r
\name{bubble_plot_by_haul_indexes}
\alias{bubble_plot_by_haul_indexes}
\title{Bubble plot of abundance and biomass indices by haul}
\usage{
bubble_plot_by_haul_indexes(
  mTATB,
  map_lim,
  depth_lines,
  buffer = 0,
  res = 0.1,
  wd = NA,
  save = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{mTATB}{data frame}

\item{map_lim}{coordinates limits for the plotted map}

\item{depth_lines}{vector of three depth bathymetrical lines to be plotted}

\item{buffer}{buffer to the coordinate limits in map units}

\item{res}{resolution of the bathymetrical lines}

\item{wd}{working directory}

\item{save}{boolean. If TRUE the plot is saved in the user defined working directory (wd)}

\item{verbose}{boolean. If TRUE a message is printed}
}
\description{
The function generates bubble plot of abundance and biomass indices by haul
}
