

#' Hamming distance between two sequences.
#'
#' @param x Character vector. Sequence 1.
#' @param y Character vector. Sequence 2.
#' @param sites List of integer vectors. Positions where the differences are calculated.
#' For example,  \code{sites = list( c(1:5), c(300:900), 1234}.
#' If \code{sites=NULL}, then all sites are included in the distance calculation.
#'
#' @return Hamming distance at the selected sites.
#' @export
#'
dist_hamming <- function(x, y, sites) {


  if(0){
    # x = res$seq[[1]]
    # y = res$seq[[2]]
    # sites = list( c(1:5), c(300:900) )
  }

  ns = length(sites)

  # If the sequence do not have the same length,
  # complete the shortest one with dummy characters.
  nx = length(x)
  ny = length(y)
  if(nx != ny){
    if(nx < ny) x = c(x, rep('X', ny-nx))
    if(ny < nx) y = c(y, rep('X', nx-ny))
  }

  # Vector of differences (0=same, 1=different)
  h <- as.numeric(x!=y)

  if(is.null(sites)){
    res = sum(h)
  }

  tmpfct <- function(i, h, sites) {
    return(h[ sites[[i]] ])
  }

  if(!is.null(sites)){
    # Differences at the selected sites:
    dh  = lapply(1:ns, tmpfct, h=h, sites=sites)
    z   = sapply(dh, sum)
    res = sum(z)
  }

  return(res)
}


#' Temporary function used in parallel computing of Hamming distance.
#'
#' @param i Integer. Index used during parallel computing.
#' @param seqs List of sequences.
#' @param sites List of integer vectors. Positions where the differences are calculated.
#' @return A numeric vector of distances.
dist_hamming_unit <- function(i,seqs, sites) {

  n = length(seqs)
  a = numeric(length = n)

  if(i==n) return(a)

  # Calculate distances only
  # for the upper triangle matrix
  # (faster!)
  for(j in c((i+1):n)){   # i=4
    a[j] = dist_hamming(x = seqs[[i]], y = seqs[[j]], sites=sites)
  }
  return(a)
}



#' Distance matrix for multiple sequences.
#'
#' @param seqs List of sequences.
#' @param sites List of integer vectors. Positions where the differences are calculated.
#' For example,  \code{sites = list( c(1:5), c(300:900), 1234}.
#' If \code{sites=NULL} (Default), then all sites are included in the distance calculation.
#' @param ncores Integer. Number of cores used for parallel computation.
#' Default \code{ncores=1} (i.e., no parallel computation). To use all available cores, enter \code{ncores=0}.
#' @param dist.type String. Type of genetic distance. Default \code{dist.type='hamming'}.
#' Options are \code{'hamming', 'blosum'}.
#'
#' @return A distance matrix.
#' @export
#'
dist_matrix <- function(seqs, sites = NULL,
                        ncores = 1,
                        dist.type = 'hamming') {

  t1 = as.numeric(Sys.time())
  # Number of sequences
  n <- length(seqs)

  message(paste('Calculating pairwise',dist.type,'distance for',n,'sequences...'))

  # Set up parallel computation:
  NC = ifelse(ncores==0, parallel::detectCores(), ncores)
  snowfall::sfInit(parallel = NC>1, cpus = NC)
  snowfall::sfExportAll()

  if(dist.type == 'hamming')
    therows <- snowfall::sfLapply(
      x     = 1:n,
      fun   = dist_hamming_unit,
      seqs  = seqs,
      sites = sites
    )

  if(dist.type == 'blosum')
    therows <- snowfall::sfLapply(x = 1:n,
                        fun = dist_blosum_unit,
                        seqs = seqs)

  snowfall::sfStop()

  m <- matrix(data = unlist(therows),
              nrow = n, ncol=n, byrow = TRUE)
  # make the matrix symetric
  # (input for `cmdscale()` must be symetric)
  m = m + t(m)

  t2 = as.numeric(Sys.time())
  dt = round((t2-t1)/60, 2)
  msg.dt = paste0('Computing time dist_matrix: ',dt,' minutes')
  message(msg.dt)
  return(m)
}




