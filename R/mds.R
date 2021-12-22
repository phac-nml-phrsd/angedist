#' Multi-dimensional scaling
#'
#' @param m Distance matrix.
#' @param dim.mds Dimension to project to.
#'
#' @return Dataframe with MDS coordinates and associated meta-variables as additional columns.
#' @export
#'
#'

mds <- function(m, dim.mds, metavars = NULL, display.GOF = TRUE) {

  if(0){
    dim.mds = 2
    metavars = list(sname = x1$strain.name, datec = x1$date.collection)
  }

  t1 = as.numeric(Sys.time())

  a = cmdscale(d = m, k = dim.mds, list. = TRUE)

  if(display.GOF){
    message('MDS goodness-of-fit by dimensions:')
    message(paste(paste0('MDSdim = ',1:dim.mds,' : ',round(a$GOF,5),'\n'),
                  collapse = ''))
  }

  mds.coord = data.frame(a$points)

  if(!is.null(metavars)){

    # Check consistency
    chk.length = sapply(metavars, length)
    stopifnot(all(chk.length == chk.length[1]))

    # Add the meta-variables
    mds.coord = cbind(mds.coord, metavars)
  }

  res = list(
    dim.mds = dim.mds,
    metavars.name = ifelse(is.null(metavars),NULL,names(metavars)),
    gof = a$GOF,
    df = mds.coord
  )

  t2 = as.numeric(Sys.time())
  dt = (t1-t2)/60
  msg.t = paste0('MDS compute time: ',dt, 'm')
  message(msg.t)
  return(res)
}


#' Plot a MDS object
#'
#' @param mdsobj MDS object as returned by the function \code{angedist::mds()}.
#' @param color_varname String. Name of the variable used for the color aesthetic.
#' For no color, \code{color_varname = NULL} (default).
#'
#' @return a \code{ggplot2} object.
#' @export
#'
#' @import ggplot2
#'
plot_mds <- function(mdsobj, color_varname = NULL) {

  # color_varname = 'date'

  df = mdsobj$df

  g = NULL

  th = theme_bw() +
    theme(axis.ticks = element_blank(),
          axis.text = element_text(color = 'grey70'),
          panel.grid.major = element_line(color='grey95'),
          panel.grid.minor = element_line(color='grey97'),
    )

  pt.sz = 3
  pt.alpha = .7

  aes12 = aes_string(x='X1', y='X2', color = color_varname)
  aes13 = aes_string(x='X1', y='X3', color = color_varname)
  aes23 = aes_string(x='X2', y='X3', color = color_varname)

  gp = geom_point(size = pt.sz, alpha = pt.alpha)

  lbs = labs(
    # x = 'antigenic dimension 1',
    # y = 'antigenic dimension 2',
    title = 'MDS'
  )

  if(mdsobj$dim.mds == 2){
    g = ggplot(df, aes12) + gp + th + lbs
    # g
  }
  if(mdsobj$dim.mds == 3){
    tmp = list(
    g12 = ggplot(df, aes12) + gp + th + lbs,
    g13 = ggplot(df, aes13) + gp + th + lbs,
    g23 = ggplot(df, aes23) + gp + th + lbs
    )
    g = patchwork::wrap_plots(tmp, guides = 'collect')
  }

  return(g)
}
