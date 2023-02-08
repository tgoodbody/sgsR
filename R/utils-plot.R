#' Plot
#'
#' @inheritParams strat_breaks
#' @inheritParams extract_metrics
#' @param dfc data.frame. Values for mraster and mraster2
#' @param coordsgrps List. Cartesian coordinates of each strata
#' @family plot
#' @name plot
#' @return Scatter plot of available raster cells coloured and delineated by stratum.
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
                      samp = 0.01) {
  #--- sample cells based on 'samp' parameter ---#

  dfc <- dfc %>%
    stats::na.omit() %>%
    dplyr::group_by(class) %>%
    dplyr::slice_sample(prop = samp)

  #--- plot individual cells coloured by associated class ---#

  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = {{ mraster }}, y = {{ mraster2 }})) +
    ggplot2::geom_point(data = dfc, alpha = 0.3, ggplot2::aes(color = as.factor(class))) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")

  #--- add class boundary boxes to delineate class extents ---#

  for (i in 1:nrow(coordsgrps)) {
    data <- coordsgrps$data[[i]]

    p <- p + ggplot2::geom_rect(
      data = data,
      ggplot2::aes(
        xmin = min({{ mraster }}),
        xmax = max({{ mraster }}),
        ymin = min({{ mraster2 }}),
        ymax = max({{ mraster2 }})
      ),
      colour = "black",
      fill = NA
    )
  }

  return(p)
}


#' Scatter Plot
#' @family plot
#' @rdname plot
#' @keywords internal
#' @param samp Numeric. Determines proportion of cells to plot
#' @param reverse Logical. Reverse x and y axis
#' @note Population in viridis and samples in red.
#' @export

plot_scatter <- function(mraster,
                         existing,
                         reverse = FALSE,
                         samp = 0.01) {
  .data <- NULL

  #--- Error management ---#
  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster.", call. = FALSE)
  }

  if (length(names(mraster)) == 1) {
    stop("Only 1 layer in `mraster` when 2 are needed.", call. = FALSE)
  } else if (length(names(mraster)) > 2) {
    message("More than 2 layers in `mraster`. Only first 2 layers will be used.")
  }

  if (!inherits(existing, "sf")) {
    stop("'existing' must be an sf object.", call. = FALSE)
  }

  if (!is.logical(reverse)) {
    stop("'reverse' must be type logical.", call. = FALSE)
  }

  if (!is.logical(reverse)) {
    stop("'reverse' must be type logical.", call. = FALSE)
  }

  if (!is.numeric(samp)) {
    stop("'samp' must be type numeric.", call. = FALSE)
  }

  if (samp > 1 | samp < 0) {
    stop("'samp' must be > 0 <= 1.", call. = FALSE)
  }

  #--- extract values to sample ---#

  samples <- extract_metrics(mraster = mraster, existing = existing)

  #--- extract raster values ---#

  vals <- terra::as.data.frame(mraster[[1:2]], na.rm = TRUE) %>%
    dplyr::slice_sample(prop = samp)

  #--- plot individual cells coloured by associated class ---#

  x <- ggplot2::sym(names(vals[1]))
  y <- ggplot2::sym(names(vals[2]))

  if (isFALSE(reverse)) {
    p <- ggplot2::ggplot(data = vals, mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
      ggplot2::geom_bin2d(bins = 70) +
      ggplot2::scale_fill_continuous(type = "viridis") +
      ggplot2::geom_point(data = samples, color = "red") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none")
  } else {
    p <- ggplot2::ggplot(data = vals, mapping = ggplot2::aes(x = .data[[y]], y = .data[[x]])) +
      ggplot2::geom_bin2d(bins = 70) +
      ggplot2::scale_fill_continuous(type = "viridis") +
      ggplot2::geom_point(data = samples, color = "red") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none")
  }

  #--- return plot ---#

  return(p)
}
