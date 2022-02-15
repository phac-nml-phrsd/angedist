
#' @title Clean strain name.
#' @description  Remove underscores and spaces from location of strain names
#'    A/Hong Kong/... ==> A/HongKong/...
#'    A/Hong_Kong/... ==> A/HongKong/...
#'
#' @param x String. The unprocessed strain name.
#'
#' @return String. Processed strain name.
#'
clean_strain_names_location_seq <- function(x) {
  # only one space:
  y1 = gsub('([AB]/\\w+)\\s(\\w+/)', '\\1\\2', x)
  # two spaces:
  y2 = gsub('([AB]/\\w+)\\s(\\w+)\\s(\\w+/)', '\\1\\2\\3', y1)
  # one underscore:
  y3 = gsub('([AB]/\\w+)_(\\w+/)', '\\1\\2', y2)
  # two underscores:
  y4 = gsub('([AB]/\\w+)_(\\w+)_(\\w+/)', '\\1\\2\\3', y3)

  return(y4)
}


#' @title Regex pathogen code in FASTA header.
#' @description Return the pathogen code as written in the FASTA header.
#'
#' @param p String. Pathogen name.
#' @param seq.source String. Source of sequence data (e.g., GISAID)
#'
#' @return String for a REGEX search.
pathogen_code <- function(p, seq.source) {

  if(0){
    p = 'influenza'
    seq.source = 'GISAID'
  }

  res = NULL

  if(seq.source == 'GISAID'){
    if(p == 'influenza') res = '[AB]'
    if(grepl('SARSCOV2',p))  res = 'hCoV-19'
    if(p == 'SARSCOV2spike')  res = 'Spike|hCoV-19'
  }
  return(res)
}

get_strainname <- function(x, pathogen, seq.source) {
  if(0){
    # pathogen = 'SARSCOV2spike'
    # pathogen = 'influenza'
    # seq.source = 'GISAID'
    # x = h[3]
    # x
  }
  pc  = pathogen_code(pathogen, seq.source)

  # positions of the vertical bars
  b = unlist(gregexpr('|',text = x, fixed = TRUE))

  i = ifelse(pathogen =='SARSCOV2spike',2,1)
  k = ifelse(pathogen =='SARSCOV2spike',8,2)
  z = substr(x, start = k, stop = b[i])
  z = stringr::str_remove_all(z, '\\|')

  res = unname(trimws(z))
  res = clean_strain_names_location_seq(res)

  return(res)
}


extract_collection_date <- function(s) {

  # --- Exact date is known ---

  d = stringr::str_extract(
    string = s,
    pattern = "\\|\\d{0,2}\\/?\\d{2}\\/\\d{4}\\||\\|\\s*\\d{4}-\\d{2}-\\d{2}\\s*\\|")

  dd = gsub(pattern = '|', replacement = '', x = d, fixed=TRUE)
  dd =  trimws(dd)

  # --- Unknown exact date ---

  # For the unknown date, GISAID provide the year
  idx.u = grepl('unknown', s)
  u = s[idx.u]

  # Retrieve the year
  tmp = stringr::str_extract(u, '\\|\\s*\\d{4}\\s+\\(Month and day unknown\\)\\s+\\|')
  tmp = stringr::str_extract(tmp, '\\d{4}')

  # Choose middle of the year, arbitrarily
  du = paste0(tmp,'-07-01')

  # Overwrite when data was unknown:
  dd[idx.u] <- du

  return(dd)
}

#' Get the collection date from a FASTA header.
#'
#' @param h A character vector of FASTA headers.
#'
#' @return A list of collection dates.
#'
get_collection_date <- function(h) {

  # Dates in character format:
  dc = extract_collection_date(h)

  # Check the date format by checking
  # the first two left-most digits.
  ff    = as.numeric( stringr::str_extract(dc, '^\\d{2}') )
  maxff = max(ff, na.rm = TRUE)

  date.format = NA
  if(maxff==12) date.format = 'mdy'
  if(maxff==31) date.format = 'dmy'
  if(maxff==20) date.format = 'ymd'

  if(is.na(date.format)){
    warning('Unable to determine the date format of sequence collection date!')
    d = rep(NA, length(h))
  }

  # Convert to Date objects:
  if(date.format=='mdy') d = lubridate::mdy(dc)
  if(date.format=='dmy') d = lubridate::dmy(dc)
  if(date.format=='ymd') d = lubridate::ymd(dc)

  # Force the output to be in the YMD format:
  res = lubridate::ymd(d)

  return(res)
}



get_location <- function(x, pathogen, seq.source) {

  if(0){
    # x = sname[[1]]
    # pathogen = 'influenza'
    # pathogen = 'SARSCOV2spike'
    # seq.source = 'GISAID'
  }
  pc  = pathogen_code(pathogen, seq.source)

  if(pathogen == 'SARSCOV2spike') {
    pc = stringr::str_remove(pc, '^Spike\\\\\\|') # TODO: change this
  }
  ss  = paste0(pc, '/\\w+')
  tmp = stringr::str_extract(x, ss)
  res = stringr::str_remove(tmp,  paste0(pc, '/'))

  return(unname(as.character(res)))
}


#' Retrieve the accession/ID number from the FASTA header.
#'
#' @param x  String. FASTA header.
#' @param seq.source String. Source of sequencing data (e.g., GISAID).
#'
#' @return String for the accession/ID number
#'
get_accession_number <- function(x, seq.source) {
  a = ''
  if(seq.source == 'GISAID') a = 'EPI_ISL_'

  r = stringr::str_extract(x, paste0(a, '\\d+'))
  return(r)
}


#' Read and import genetic sequences from files.
#'
#' @param path String. Path to the file storing the genetic sequences.
#' @param prms List. Parameters for importation. Must specify \code{seq.source},
#'  \code{pathogen}, \code{seq.type} (\code{AA} or \code{DNA}).
#'
#'
#' @return A \code{SeqFastaXXX} object as defined by the package \code{seqinr}.
#' @export
#'
#' @importFrom seqinr read.fasta
#'
import_seqs <- function(path, prms) {

  t1 = as.numeric(Sys.time())

  message(paste0('Reading genetic sequence at: ', path, ' ...'),
          appendLF = FALSE)

  # Retrieve the sequence objects from the FASTA file
  s = seqinr::read.fasta(file = path,
                         seqtype = prms$seq.type)
  s = unname(s)

  # Retrieve FASTA sequence header
  h = sapply(s, attr, 'Annot')

  # Full strain/virus name
  sname = unname(sapply(h, get_strainname,
                        pathogen = prms$pathogen,
                        seq.source = prms$seq.source))

  # Collection dates
  dc = get_collection_date(h)

  # Collection location
  loc = sapply(sname, get_location,
               pathogen   = prms$pathogen,
               seq.source = prms$seq.source,
               USE.NAMES  = FALSE)

  # Accession number
  an = sapply(h, get_accession_number,
              seq.source = prms$seq.source, USE.NAMES = FALSE)

  res = list(
    seq = s,
    strain.name = sname,
    date.collection = dc,
    location = loc,
    accession.num = an,
    header = h,
    prms = prms
  )

  t2 = as.numeric(Sys.time())
  dt = round(t2-t1,1)
  msg.t = paste0(' done (time = ',dt,' sec).')
  message(msg.t)
  return(res)
}
