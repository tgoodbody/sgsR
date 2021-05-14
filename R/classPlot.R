#' Scatter plot of strata using ggplot
#'
#' @inheritParams strat_quantiles
#' @param dfc data.frame. Values for metric1 and metric2
#' @param coordsgrps List. Cartesian coordinates of each strata
#' @param samp Numeric. Determines proportion of cells to plot 
#' for strata visualization. Lower values reduce processing time.
#' 
#' @importFrom magrittr %>%
#' @importFrom methods is
#'


classPlot <- function(dfc,
                      coordsgrps,
                      metric,
                      metric2,
                      samp){

  metric <- ggplot2::sym(metric)
  metric2 <- ggplot2::sym(metric2)
  
  #--- sample cells based on 'samp' parameter ---#
  
  dfc <- dfc %>%
    stats::na.omit() %>%
    dplyr::group_by(class) %>%
    dplyr::slice_sample(prop = samp)
  
  #--- plot individual cells coloured by associated class ---#

  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = !!metric, y = !!metric2))+
    ggplot2::geom_point(data = dfc, alpha = 0.3, ggplot2::aes(color = as.factor(class)))+
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
  
  #--- add class boundary boxes to delineate class extents ---#

  for(i in 1:nrow(coordsgrps)){
    
    data <- coordsgrps$data[[i]]
    
    p <- p + ggplot2::geom_rect(data = data,
                      ggplot2::aes(xmin = min(!!metric),
                           xmax = max(!!metric),
                           ymin = min(!!metric2),
                           ymax = max(!!metric2)),
                       colour = "black",
                       fill = NA)
  }

  return(p)
}
