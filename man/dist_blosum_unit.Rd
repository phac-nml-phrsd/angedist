% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dist.R
\name{dist_blosum_unit}
\alias{dist_blosum_unit}
\title{Helper function for parallel computing.}
\usage{
dist_blosum_unit(i, seqs, sites, logistic.scale = 15)
}
\arguments{
\item{i}{Integer. Iteration number.}

\item{seqs}{List of sequences.}

\item{sites}{List of integer vectors. Positions where the differences are calculated.
For example,  \code{sites = list( c(1:5), c(300:900), 1234}.
Default \code{logistic.scale = 15}.}

\item{logistic.scale}{Numeric. Scale parameter for the
logistic function when converting similarity into a distance.}
}
\value{
Numeric. BLOSUM distance.
}
\description{
Helper function for parallel computing.
}
