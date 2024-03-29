% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_dictionary.r
\name{check_dictionary}
\alias{check_dictionary}
\title{Check dictionary (RoME)}
\usage{
check_dictionary(
  ResultData,
  Field,
  Values,
  year,
  wd = NA,
  suffix,
  verbose = FALSE
)
}
\arguments{
\item{ResultData}{data frame in MEDITS tables}

\item{Field}{field of the table to be checked}

\item{Values}{vector of the allowed values}

\item{year}{reference year for the analysis}

\item{wd}{working directory}

\item{suffix}{name of the log file}

\item{verbose}{boolean. If TRUE messages are promted in the console}
}
\description{
The function checks whether the values contained in specific fields are consistent with the allowed values of the dictionaries.
}
