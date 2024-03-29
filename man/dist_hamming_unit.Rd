% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dist.R
\name{dist_hamming_unit}
\alias{dist_hamming_unit}
\title{Temporary function used in parallel computing of Hamming distance.}
\usage{
dist_hamming_unit(i, seqs, sites)
}
\arguments{
\item{i}{Integer. Index used during parallel computing.}

\item{seqs}{List of sequences.}

\item{sites}{List of integer vectors. Positions where the differences are calculated.}
}
\value{
A numeric vector of distances.
}
\description{
Temporary function used in parallel computing of Hamming distance.
}
