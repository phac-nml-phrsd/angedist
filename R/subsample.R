
#' Sub-sample sequences based on their collection date
#'
#' @param sobj List of sequences as returned by \code{import_seqs()}.
#' @param n.subsample Integer. Size of the subsample.
#' @param years.exception Integer vector. Exception years
#' (sequences from those years will not be sub-sampled).
#' If no exception, enter \code{NULL}.
#' @param verbose Logical.
#'
#' @return A list of sequences that have been sub-sampled.
#'
#' @export
#'
#' @import lubridate
#'
subsample_date <- function(sobj,
                           n.subsample,
                           years.exception = NULL,
                           verbose = FALSE) {

  N = length(sobj$date.collection)
  y = lubridate::year(sobj$date.collection)

  if(!is.null(years.exception)){
    # indices to keep because their collection year is an exception:
    idx.xcpt = which(y %in% years.exception)
    # indices to sample from:
    idx.smpl = c(1:N)[-idx.xcpt]
  }
  if(is.null(years.exception)){
    idx.xcpt = NULL
    idx.smpl = c(1:N)
  }

  idx = sample(idx.smpl, size = n.subsample)
  # merge the subsampled indices and
  # the ones satisfying the exception:
  idx.subsmpl = c(idx, idx.xcpt)

  # build the subsampled sequence object:
  res = list(
    seq             = sobj$seq[idx.subsmpl],
    strain.name     = sobj$strain.name[idx.subsmpl],
    date.collection = sobj$date.collection[idx.subsmpl],
    location        = sobj$location[idx.subsmpl],
    accession.num   = sobj$accession.num[idx.subsmpl]
  )
  res = restore_unchanged(res, sobj)

  if(verbose){
    msg = paste0('Subsampling among ',N, ' sequences:\n',
                 'Subsample size requested = ', n.subsample,'\n',
                 'Addition from exceptions = ',length(idx.xcpt),'\n',
                 'Final subsample size     = ', n.subsample + length(idx.xcpt))
    message(msg)
  }
  return(res)
}

