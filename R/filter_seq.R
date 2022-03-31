

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
#' @param except.strain.name String vector of strain name to exclude from the
#' filtering. That means the sequences returned may not be all from the
#' specified locations (if the exceptions are from other locations).
#'
#' @return A list of sequence objects.
#' @export
#'
filter_location <- function(sobj, loc, except.strain.name=NULL) {

  if(0){  # --- DEBUG
    loc = c('Ontario', 'Manitoba')
    except.strain.name = c('A/Nagasaki/97/95',
                           'A/Denmark/41/03')
  }

  p = tolower(loc)
  q = tolower(sobj$location)

  # string pattern for `grepl()`:
  pp = paste(p, collapse = '|')

  # return the indices where the locations were found:
  idx = which(grepl(pp, q))

  # No locations found
  if(length(idx)==0) res = NULL


  # Add sequences from the exception
  if(!is.null(except.strain.name)){
    xsn = tolower(except.strain.name)
    sn  = tolower(sobj$strain.name)
    idx.xsn = which(sn %in% xsn)

    if(length(idx.xsn)==0){
      warning(paste('None of the',length(except.strain.name),
                    'strain name exceptions were found in `filter_location()`.'))
    }
    if(length(idx.xsn) < length(except.strain.name) & (length(idx.xsn)>0) ){
      warning(paste('Only',length(idx.xsn),'of the',
                    length(except.strain.name),
                    'strain name exceptions were found in `filter_location()`.'))
    }

    # add to the filtered locations:
    idx = c(idx, idx.xsn)
  }
  if(length(idx)>0)  res = filter_index(sobj, idx)

  return(res)
}
