% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mds.R
\name{mds}
\alias{mds}
\title{Multi-dimensional scaling}
\usage{
mds(m, dim.mds, metavars = NULL, display.GOF = TRUE)
}
\arguments{
\item{m}{Distance matrix as returned from the function \code{dist_matrix()}.}

\item{dim.mds}{Dimension to project to.}
}
\value{
Dataframe with MDS coordinates and associated meta-variables as additional columns.
}
\description{
Multi-dimensional scaling
}
