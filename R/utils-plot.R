#' Plot
#'
#' @inheritParams strat_breaks
#' @param dfc data.frame. Values for mraster and mraster2
#' @param coordsgrps List. Cartesian coordinates of each strata
#' @param samp Numeric. Determines proportion of cells to plot
#' @family plot
#' @name plot
NULL

#' Class Plot
#' @family plot
#' @rdname plot
#' @keywords internal
#' @export

classPlot <- function(dfc,
                      coordsgrps,
                      mraster,
                      mraster2,
                      samp) {

  #--- if the variables are not symbols convert them ---#

  if (is.character(mraster)) {
    mraster <- ggplot2::ensym(mraster)
  }

  if (is.character(mraster2)) {
    mraster2 <- ggplot2::ensym(mraster2)
  }

  #--- sample cells based on 'samp' parameter ---#

  dfc <- dfc %>%
    stats::na.omit() %>%
    dplyr::group_by(class) %>%
    dplyr::slice_sample(prop = samp)

  #--- plot individual cells coloured by associated class ---#

  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = !!mraster, y = !!mraster2)) +
    ggplot2::geom_point(data = dfc, alpha = 0.3, ggplot2::aes(color = as.factor(class))) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")

  #--- add class boundary boxes to delineate class extents ---#

  for (i in 1:nrow(coordsgrps)) {
    data <- coordsgrps$data[[i]]

    p <- p + ggplot2::geom_rect(
      data = data,
      ggplot2::aes(
        xmin = min(!!mraster),
        xmax = max(!!mraster),
        ymin = min(!!mraster2),
        ymax = max(!!mraster2)
      ),
      colour = "black",
      fill = NA
    )
  }

  return(p)
}
