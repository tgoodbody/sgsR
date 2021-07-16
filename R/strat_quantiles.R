#' Quantiles stratification
#'
#' @description Stratify metric raster using metric quantiles.
#' @family stratify functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams strat_breaks
#' @param metric Character. Name of primary metric to stratify. If
#' \code{mraster} is has 1 layer it is taken as default.
#' @param metric2 Character. Name of secondary metric to stratify.
#' @param nStrata2 Numeric.  Number of secondary strata within \code{nStrata}.
#' @param samp Numeric. Determines proportion of cells to plot in scatterplot (see \code{values})
#' for strata visualization. Lower values reduce visualization time.
#'
#' @importFrom methods is
#'
#' @return Returns an output stratification \code{spatRaster} or a list when \code{details = TRUE}.
#'  
#' When a list is returned:
#' \enumerate{
#' \item \code{details} is a list output of the \code{\link[stats]{prcomp}} function
#' \item \code{raster} is a stratified \code{spatRaster} based on quantiles
#' \item \code{scatter} is a \code{ggplot} histogram / scatter plot object (depends on whether metric2 was supplied).
#' Histogram shows distribution and break points while scatter plot shows colour coded and strata boundaries.
#' }
#'
#' @export

strat_quantiles <- function(mraster,
                            metric = NULL,
                            metric2 = NULL,
                            nStrata,
                            nStrata2 = NULL,
                            plot = FALSE,
                            details = FALSE,
                            samp = 1,
                            filename = NULL,
                            overwrite = FALSE,
                            ...) {

  #--- Set global vars ---#

  class1 <- class2 <- NULL

  #--- error handling ---#
  if (!inherits(mraster, "SpatRaster")) {
    stop("all specified bands must be type SpatRaster", call. = FALSE)
  }

  if (!is.numeric(nStrata)) {
    stop("'nStrata' must be type numeric")
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical")
  }

  if (!is.numeric(samp)) {
    stop("'samp' must be type numeric")
  }

  if (!is.logical(details)) {
    stop("'details' must be type logical")
  }

  if (is.null(metric2)) {
    if (!is.null(nStrata2)) {
      message("You are stratifying with only 1 metric but specified 'nStrata2' - ignoring.")
    }

    #--- if there is only 1 metric in the raster use it as default ---#

    if (terra::nlyr(mraster) == 1) {

      #--- Extract values from mraster ---#

      vals <- terra::values(mraster)

      #--- set name of raster band to 'metric' ---#

      metric <- names(mraster)
    } else {
      if (is.null(metric)) {
        stop(" multiple layers detected in 'mraster'. Please define a 'metric' to stratify")
      }

      if (!is.character(metric)) {
        stop("'metric' must be type character")
      }

      if (any(!metric %in% names(mraster))) {
        stop(paste0("'mraster' must have an attribute named ", metric))
      }

      #--- Extract values from mraster ---#

      vals <- terra::subset(mraster, metric) %>%
        terra::values()
    }

    vals[!is.finite(vals)] <- NA

    #--- Determine index of each cell so to map values correctly without NA ---#

    idx <- !is.na(vals)

    #--- Remove NA / NaN / Inf values ---#

    df <- vals %>%
      as.data.frame() %>%
      dplyr::filter(!is.na(.))

    metric <- ggplot2::ensym(metric)

    #--- Split metric distribution in to number specified by 'breaks' ---#

    dfc <- df %>%
      dplyr::mutate(class = dplyr::ntile(!!metric, nStrata))

    #--- convert back to original mraster extent ---#

    vals[idx] <- dfc$class

    #--- set newly stratified values ---#

    rout <- terra::setValues(mraster[[1]], vals)
    names(rout) <- "strata"
  }

  if (!is.null(metric2)) {
    if (!is.character(metric2)) {
      stop("'metric2' must be type character")
    }

    if (is.null(nStrata2)) {
      stop("If using 2 metrics to stratify, 'nStrata2' must be defined")
    }

    if (any(!metric2 %in% names(mraster))) {
      stop(paste0("'mraster' must have an attribute named ", metric2))
    }

    #--- Extract values from mraster ---#

    vals <- terra::subset(mraster, c(metric, metric2)) %>%
      terra::values()

    vals[!is.finite(vals)] <- NA

    #--- Determine index of each cell so to map values correctly without NA ---#

    idx <- is.finite(vals[, 1]) & is.finite(vals[, 2])

    #--- Remove NA / NaN / Inf values ---#

    df <- vals %>%
      as.data.frame() %>%
      dplyr::filter(!is.na(.))

    metric <- ggplot2::ensym(metric)
    metric2 <- ggplot2::ensym(metric2)

    #--- Split metric distribution in to number specified by 'breaks' ---#

    dfc <- df %>%
      #--- define nStrata classes ---#
      dplyr::mutate(class1 = dplyr::ntile(!!metric, nStrata)) %>%
      #--- group by class to sub stratify ---#
      dplyr::group_by(class1) %>%
      #--- define nStrata2 classes ---#
      dplyr::mutate(class2 = dplyr::ntile(!!metric2, nStrata2)) %>%
      #--- combine classes ---#
      dplyr::group_by(class1, class2) %>%
      #--- establish newly formed unique class ---#
      dplyr::mutate(class = dplyr::cur_group_id())

    #--- convert back to original mraster extent ---#

    vals[, 1][idx] <- dfc$class

    #--- set newly stratified values ---#

    rout <- terra::setValues(mraster[[1]], vals[, 1])
    names(rout) <- "strata"
  }

  if (isTRUE(plot)) {
    if (samp > 1 | samp < 0) {
      stop("'samp' must be between 0 and 1")
    }

    if (is.null(metric2)) {

      #--- output histogram with quantile breaks ---#

      #--- determine numeric breaks ---#

      breaks <- dfc %>%
        dplyr::group_by(class) %>%
        dplyr::summarize(breaks = max(!!metric)) %>%
        dplyr::select(breaks) %>%
        as.data.frame()

      breaks <- breaks[1:(nStrata - 1), ]

      df.p <- dfc %>%
        dplyr::select(metric)

      #--- plot histogram of metric with associated break lines ---#

      p <- ggplot2::ggplot(df.p, ggplot2::aes(!!metric)) +
        ggplot2::geom_histogram() +
        ggplot2::geom_vline(xintercept = breaks, linetype = "dashed") +
        ggplot2::ggtitle(paste0(metric, " histogram with defined breaks"))

      print(p)
    } else {

      #--- set up colour palette ---#

      terra::plot(rout, main = "Classes")

    }

    #--- write file to disc ---#

    if (!is.null(filename)) {
      terra::writeRaster(rout, filename, overwrite = overwrite, ...)
    }

    #--- Output based on 'details' to return raster alone or list with details ---#

    if (isTRUE(details)) {
      
      #--- create classplot summary ---#
      
      coordsgrps <- dfc %>%
        dplyr::group_by(class) %>%
        dplyr::arrange(class) %>%
        tidyr::nest() %>%
        dplyr::ungroup()
      
      p <- classPlot(
        dfc = dfc,
        coordsgrps = coordsgrps,
        metric = metric,
        metric2 = metric2,
        samp = samp
      )

      #--- output metrics details along with stratification raster ---#

      out <- list(
        details = dfc,
        raster = rout,
        plot = p
      )

      return(out)
    } else {

      #--- just output raster ---#

      return(rout)
    }
  }

  return(rout)
}
