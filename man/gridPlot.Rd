% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gridPlot.R
\name{gridPlot}
\alias{gridPlot}
\title{Plot DGP Parameters}
\usage{
gridPlot(data, plot = c("posterior", "conditional"), cutoff = 0.05)
}
\arguments{
\item{data}{A data frame generated by \code{gridSearch}.}

\item{plot}{If \code{posterior}, plots a density plot for each non-fixed parameter. If \code{conditional}, plots the accepted samples by different balance levels.}

\item{cutoff}{The KS distance cutoff for accepting the ABC samples.}
}
\value{
A \code{ggplot2} object.
}
\description{
This function plots the distribution of accepted samples.
}
