

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

  tmp_hamming <- function(i, h, sites) {
    return(h[ sites[[i]] ])
  }

  if(!is.null(sites)){
    # Differences at the selected sites:
    dh  = lapply(1:ns, tmp_hamming, h=h, sites=sites)
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



#' @title BLOSUM distance between two sequences.
#'
#' @description Calculate the genetic distance according to
#' the BLOSUM metric. Here, only the BLOSUM80 is implemented.
#' BLOSUM80 definition: https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM80
#'
#'
#' @param x Character vector representing the amino acids sequence.
#' @param y Character vector representing the amino acids sequence.
#' @param logistic.scale Numeric. Scale parameter for the logistic function
#' when converting similarity into a distance.
#' Default \code{logistic.scale = 15}.
#'
#' @return Numeric. BLOSUM80 distance between \code{x} and \code{y}.
#' @export
#'
#' @examples
#' x = c('A','P','P','L','E','S')
#' y = c('A','P','M','L','E','S')
#' dist_blosum(x,y,  sites = list(c(1:2), 4:6))
#' dist_blosum(x,y,  sites = NULL)

dist_blosum <- function(x, y, sites = NULL,
                        logistic.scale = 15) {

  n  = length(x)

  xx = x; yy=y
  if(!is.null(sites)){
    s = unlist(sites)
    xx = x[s]
    yy = y[s]
  }

  m  = blosum80[xx,yy]  # `blosum80` is defined in `blosum.R`
  similarity.score = sum(diag(m)) / n

  # the BLOSUM are log-odds similarities.
  # Use logistic function to translate in distance:
  d = 1 - plogis(q = -similarity.score, scale = logistic.scale)
  return(d)
}






#' @title Helper function for parallel computing.
#'
#' @param i Integer. Iteration number.
#' @param seqs List of sequences.
#' @param logistic.scale Numeric. Scale parameter for the
#' logistic function when converting similarity into a distance.
#' @param sites List of integer vectors. Positions where the differences are calculated.
#' For example,  \code{sites = list( c(1:5), c(300:900), 1234}.
#' Default \code{logistic.scale = 15}.
#' @return Numeric. BLOSUM distance.
#'
dist_blosum_unit <- function(i, seqs, sites,
                             logistic.scale = 15) {

  n = length(seqs)
  a = numeric(length = n)

  if(i==n) return(a)

  for(j in c((i+1):n)){   # i=4
    a[j] = dist_blosum(x = seqs[[i]], y = seqs[[j]], sites=sites )
  }
  return(a)
}


#' Distance matrix for multiple sequences.
#'
#' @param sobj List of sequences object as returned by \code{import_seqs()}.
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
dist_matrix <- function(sobj,
                        sites = NULL,
                        ncores = 1,
                        dist.type = 'hamming') {

  t1 = as.numeric(Sys.time())

  # Number of sequences
  seqs = sobj$seq
  n    = length(seqs)

  message(paste('\nCalculating pairwise',dist.type,
                'distance for',n,'sequences...'))

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
    therows <- snowfall::sfLapply(
      x     = 1:n,
      fun   = dist_blosum_unit,
      seqs  = seqs,
      sites = sites
    )

  snowfall::sfStop()

  m <- matrix(data = unlist(therows),
              nrow = n, ncol=n, byrow = TRUE)
  # make the matrix symetric
  # (input for `cmdscale()` must be symmetric)
  m = m + t(m)

  rownames(m) <- sobj$strain.name
  colnames(m) <- sobj$strain.name

  t2 = as.numeric(Sys.time())
  dt = round((t2-t1)/60, 2)
  msg.dt = paste0('Computing time dist_matrix: ',dt,' minutes')
  message(msg.dt)
  return(m)
}

#' Locate the most frequent starting position
#' of sub-sequence among multiple sequences.
#'
#' @param sobj Sequence object as returned by \code{import_seqs()}.
#' @param subseq String. Sub-sequence to locate.
#'
#' @return The most frequent location.
#' @export
#'
locate_postion <- function(sobj, subseq) {

  # translate into string type:
  s = lapply(sobj$seq, paste0, collapse='')
  # locate position for al sequences:
  a = unlist(gregexpr(subseq, s))

  n.found = sum(a>0)
  if(n.found ==0){
    return(-1)
  }

  # only keep the sequences that have the sub-sequence:
  p = a[a>0]

  # if multiple positions, take the most frequent:
  frq = as.matrix(table(p))
  imax = which.max(frq)

  res = as.numeric(rownames(frq)[imax])

  if(imax>1){
    fmax = frq[imax,1] / sum(frq[,1])
    warning(paste0('Subsequence ',subseq,
                   ' located at more than one position. ',
                   'Using most frequent (freq=',
                   round(fmax,4),
                   ') position (pos=',
                   res,').'))
  }

  return(res)
}



