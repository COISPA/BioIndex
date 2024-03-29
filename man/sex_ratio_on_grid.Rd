% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sex_ratio_on_grid.r
\name{sex_ratio_on_grid}
\alias{sex_ratio_on_grid}
\title{Plot sex ratio spatial distribution}
\usage{
sex_ratio_on_grid(
  mTATBsp,
  depth,
  wd,
  map_range,
  threshold = 30,
  verbose = FALSE,
  save = FALSE
)
}
\arguments{
\item{mTATBsp}{spatial mTATB}

\item{depth}{reference depth range}

\item{wd}{working directory}

\item{map_range}{range of coordinates for the map}

\item{threshold}{minimum number of individuals per haul}

\item{verbose}{boolean. If TRUE messages are prompted in the console}

\item{save}{boolean. If TRUE the outputs are saved in the local folder}
}
\description{
Plot sex ratio spatial distribution
}
