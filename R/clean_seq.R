

restore_unchanged <- function(new, original) {
  new$prms <- original$prms
  return(new)
}


#' @title Filter sequences based on their length
#'
#' @param x Sequence list, as returned by the function \code{angedist::import_seqs()}
#' @param length.min Integer. Minimum length (inclusive).
#' @param length.max Integer. Maximum length (inclusive).
#' @param verbose Logical. Display additional information about the filtering process.
#'
#' @return A list of sequences that all comply to the length specifications.
#' @export
#'
#'
filter_seq_length <- function(x, length.min, length.max, verbose=FALSE) {
  z    = sapply(x$seq, FUN = length)
  idx  = (length.min <= z & z <=  length.max)
  idx2 = c(1:length(idx))[idx]
  y    = lapply(x, '[', idx2)

  # Restore unchanged components:
  y = restore_unchanged(y, x)

  if(verbose) message(paste('--> sequences removed because of length:',length(z)-sum(idx)))
  return(y)
}


#' Remove low-quality sequences.
#'
#' @param x List of sequences.
#' @param maxreltol Numerical. Maximum relative tolerance.
#' @param unwanted Character vector of unwanted characters in the sequence. Default is set to \code{c('X','x','-')}.
#' @param verbose Logical. Display additional information about the filtering process.
#'
#' @return A list of sequences that all comply with the tolerance for missing/aberrant sata.
#' @export
#'
remove_lowquality_seqs <- function(x,
                                   maxreltol ,
                                   unwanted = c('X','x','-'),
                                   verbose = FALSE){

  # proportion unwanted
  p_unwanted <- function(z, unwanted) {
    return(mean(z %in% unwanted))
  }
  p = sapply(x$seq, FUN = p_unwanted, unwanted=unwanted)
  keepidx = (p <= maxreltol)
  y = lapply(x, '[', keepidx)

  # Restore unchanged components:
  y = restore_unchanged(y, x)

  if(verbose){
    message(paste('--> sequences removed because of missing/aberrant:',
                  length(x$seq)-sum(keepidx)))
    m = data.frame(name_removed = x$strain.name[!keepidx],
                   prop_unwanted = p[!keepidx])
    print(m)
  }
  return(y)
}




