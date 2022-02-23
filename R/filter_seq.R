

#' Title
#'
#' @param sobj List as returned by the function \code{import_seqs()}.
#' @param idx Integer vector representing the indices of the sequences selected.
#'
#' @return  A list of sequence objects.
#' @export
#'
filter_index <- function(sobj, idx) {

  res = list(
    seq             = sobj$seq[idx],
    strain.name     = sobj$strain.name[idx],
    date.collection = sobj$date.collection[idx],
    location        = sobj$location[idx],
    accession.num   = sobj$accession.num[idx],
    header          = sobj$header[idx],
    prms            = sobj$prms
  )
  return(res)
}

#' Filter sequences based on the location name.
#'
#' @param sobj List as returned by the function \code{import_seqs()}.
#' @param loc String vector of pattern to be searched for the locations.
#' If more than one location, will return all sequences that match any
#' of the the location (logical "OR")
#'
#' @return A list of sequence objects.
#' @export
#'
filter_location <- function(sobj, loc) {

  # loc = c('Ontario', 'Manitoba')
  # loc = c('qwpqwpqwp')

  p = tolower(loc)
  q = tolower(sobj$location)

  # string pattern for `grepl()`:
  pp = paste(p, collapse = '|')

  # return the indices where the locations were found:
  idx = which(grepl(pp, q))

  if(length(idx)==0) res = NULL
  if(length(idx)>0)  res = filter_index(sobj, idx)

  return(res)
}
