#' Multi-dimensional scaling
#'
#' @param m Distance matrix as returned from the function \code{dist_matrix()}.
#' @param dim.mds Dimension to project to.
#'
#' @return Dataframe with MDS coordinates and associated meta-variables as additional columns.
#' @export
#'
#'

mds <- function(m, dim.mds, metavars = NULL, display.GOF = TRUE) {

  if(0){ # -- DEBUG
    dim.mds = 2
    metavars = NULL
    metavars = list(sname = seqs.aligned$strain.name,
                    datec = seqs.aligned$date.collection)
  }

  t1 = as.numeric(Sys.time())

  a = cmdscale(d = m, k = dim.mds, list. = TRUE)

  if(display.GOF){
    message(paste('MDS goodness-of-fit:',
                  paste(round(a$GOF,5), collapse = ';')))
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
    metavars.name = ifelse(is.null(metavars),NA,names(metavars)),
    gof = a$GOF,
    df = mds.coord
  )

  t2 = as.numeric(Sys.time())
  dt = round( (t2-t1)/60, 1)
  msg.t = paste0('MDS compute time: ',dt, ' min')
  print(msg.t)
  return(res)
}


#' Plot a MDS object
#'
#' @param mdsobj MDS object as returned by the function \code{angedist::mds()}.
#' @param color_varname String. Name of the variable used for the color aesthetic.
#' For no color, \code{color_varname = NULL} (default).
#' @param highlight Logical. Highlight MDS points that are identified as being
#' highlighted by defining them through a variable called \code{highlight}
#' in the dataframe \code{mdsobj$df}.
#' @param highlight.label Logical. Label highlighted MDS points with text defined
#' by the variable \code{highlight.label} in the dataframe \code{mdsobj$df}.
#'
#' @return A \code{ggplot2} object.
#' @export
#'
#' @import ggplot2
#' @import ggrepel
#'
plot_mds <- function(mdsobj,
                     color_varname = NULL,
                     highlight = FALSE,
                     highlight.label = FALSE) {

  # highlight = 1

  df = mdsobj$df

  if(highlight & (!'highlight' %in% names(df))){
    warning('Warning: Cannot highlight MDS points because no `highlight` variable has been defined')
    highlight = FALSE
  }
  df.highlight = filter(df, highlight==TRUE)

  g = NULL

  th = theme_bw() +
    theme(axis.ticks = element_blank(),
          axis.text = element_text(color = 'grey70'),
          panel.grid.major = element_line(color='grey95'),
          panel.grid.minor = element_line(color='grey97'),
    )

  pt.sz = 3
  pt.alpha = .7
  col.highlight = 'yellow'

  aes12 = aes_string(x='X1', y='X2', color = color_varname)
  aes13 = aes_string(x='X1', y='X3', color = color_varname)
  aes23 = aes_string(x='X2', y='X3', color = color_varname)

  gp = geom_point(size = pt.sz, alpha = pt.alpha)

  gp.h    = NULL
  gp.hlab = NULL
  if(highlight){
    gp.h = geom_point(data = df.highlight,
                              size = pt.sz+1,
                              shape = 21,
                              fill = col.highlight)
  }

  if(highlight.label){
    gp.hlab = ggrepel::geom_label_repel(data = df.highlight,
                                aes(label = highlight.label))
  }

  lbs = labs(
    # x = 'antigenic dimension 1',
    # y = 'antigenic dimension 2',
    # title = 'MDS'
  )

  if(mdsobj$dim.mds == 2){
    g = ggplot(df, aes12) + gp + gp.h + gp.hlab + th + lbs
    # g
  }
  if(mdsobj$dim.mds == 3){
    tmp = list(
    g12 = ggplot(df, aes12) + gp + gp.h + gp.hlab + th + lbs,
    g13 = ggplot(df, aes13) + gp + gp.h + gp.hlab + th + lbs,
    g23 = ggplot(df, aes23) + gp + gp.h + gp.hlab + th + lbs
    )
    g = patchwork::wrap_plots(tmp, guides = 'collect')
  }

  return(g)
}
