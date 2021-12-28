

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

  if(verbose) {
    message(paste('--> sequences removed because of length:',
                  length(z)-sum(idx)))
    print(x$strain.name[!idx])
    }

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


#' Filter sequences that have defined starting proteins.
#'
#' @param x List of sequences.
#' @param subseq String. Sequence of characters that must start the sequence.
#' For example \code{subseq='MK'} (not \code{c('M','K')}).
#' @param verbose Logical. Display additional information about
#' the filtering process.
#' @param position String. Position of the subsequence:
#' \code{"start"} or \code{"end"}.
#'
#' @return A list of sequences that all comply with the tolerance for missing/aberrant sata.
#' @export
#'
filter_subseq<- function(x, subseq,
                         position,
                         verbose = FALSE) {

  tmpfct <- function(z, subseq, position){
     # z = x$seq[[1]]
    nss = nchar(subseq)
    if(position == 'start') idx = 1:nss
    if(position == 'end') idx = ( (length(z)-nss+1):length(z) )

    ssvec = unlist(strsplit(subseq, split=''))
    chk = (z[idx] == ssvec)
    return(all(chk))
  }
  keepidx <- sapply(x$seq, FUN=tmpfct, subseq=subseq, position=position)

  y = lapply(x, '[', keepidx)
  # Restore unchanged components:
  y = restore_unchanged(y, x)

  if(verbose){
    message(paste('--> sequences removed because no',subseq,'start:',
                  length(x$seq)-sum(keepidx)))
    print(x$strain.name[!keepidx])
  }
  return(y)
}

#' @title Remove duplicated sequences.
#'
#' @param x List of sequences.
#' @param verbose Logical. Display additional information about the filtering process.
#'
#' @return List of non-duplicated sequences.
#' @export
#'
remove_duplicated <- function(x, verbose = FALSE) {

  # x = res

  v = x$strain.name
  keepidx = !duplicated(v)

  y = lapply(x, '[', keepidx)
  # Restore unchanged components:
  y = restore_unchanged(y, x)

  # TODO: optional check that the sequences are indeed the same.
  # (we only test on strain name at the moment)
  # ii = which(idx.rm)
  # v[ii]
  # v[c( (ii-3) : (ii+3) )]
  # d = which(v=="hCoV-19/USA/CT-Yale-011/2020")
  # identical(x$seq[[d[1]]], x$seq[[d[2]]])
  # mean(x$seq[[d[1]]] != x$seq[[d[2]]])

  if(verbose){
    message(paste('--> sequences removed because duplicated:',
                  sum(!keepidx) ))
    print(x$strain.name[!keepidx])
  }
  return(y)
}




#' @title Clean influenza A sequences
#'
#' @param seqs.obj List of sequences as returned from \code{import_seqs()}
#' @param maxreltol Numeric. Maximum tolerance for the
#' relative proportion of missing/aberrant proteins in any given sequence.
#' Default = 0.005.
#'
#' @return A list of sequences cleaned according to influenza A specifications.
#' @export
#'
clean_seq_influenza_A <- function(seqs.obj, maxreltol=0.005) {
  res = seqs.obj %>%
    filter_seq_length(565, 566) %>%
    remove_lowquality_seqs(maxreltol = maxreltol) %>%
    filter_subseq('M', position = 'start') %>%
    filter_subseq('CI',position = 'end') %>%
    remove_duplicated()
    return(res)
}


#' Clean influenza B sequences
#'
#' @param seqs.obj List of sequences as returned from \code{import_seqs()}
#' @param maxreltol Numeric. Maximum tolerance for the
#' relative proportion of missing/aberrant proteins in any given sequence.
#' Default = 0.005.
#'
#' @return A list of sequences cleaned according to influenza A specifications.
#' @export
#'
clean_seq_influenza_B <- function(seqs.obj, maxreltol = 0.005) {
  res = seqs.obj %>%
    filter_seq_length(582, 587) %>%
    remove_lowquality_seqs(maxreltol = maxreltol) %>%
    filter_subseq('M', position = 'start') %>%
    filter_subseq('CL',position = 'end') %>%
    remove_duplicated()
  return(res)
}
