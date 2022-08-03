#' Quantiles stratification
#'
#' @description Stratify metric raster using metric quantiles.
#' @family stratify functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams strat_breaks
#' @param nStrata Numeric. Number of quantiles to stratify \code{mraster}.
#' @param nStrata2 Numeric. Number of quantiles to stratify \code{mraster2}.
#' @param samp Numeric. Determines proportion of cells to plot in scatterplot (default = \code{0.2}).
#' Lower values reduce visualization time.
#'
#' @return Returns an output stratification \code{spatRaster} or a list when \code{details = TRUE}.
#'
#' When a list is returned:
#' \enumerate{
#' \item \code{details} generates is a strata lookUp table for the stratification.
#' \item \code{raster} is a stratified \code{spatRaster} based on quantiles
#' \item \code{plot} is a \code{ggplot} histogram / scatter plot object (depends on whether metric2 was supplied).
#' Histogram shows distribution and break points while scatter plot shows colour coded and strata boundaries.
#' }
#'
#' @examples
#' #--- Load raster and existing plots---#
#' r <- system.file("extdata", "mraster.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' strat_quantiles(
#'   mraster = mr$zq90,
#'   nStrata = 4,
#'   plot = TRUE
#' )
#'
#' strat_quantiles(
#'   mraster = mr$zq90,
#'   mraster2 = mr$zsd,
#'   nStrata = 3,
#'   nStrata2 = 4
#' )
#' @author Tristan R.H. Goodbody
#'
#' @export

strat_quantiles <- function(mraster,
                            mraster2 = NULL,
                            nStrata,
                            nStrata2 = NULL,
                            plot = FALSE,
                            details = FALSE,
                            samp = 0.2,
                            filename = NULL,
                            overwrite = FALSE) {

  #--- Set global vars ---#

  class1 <- class2 <- NULL

  #--- error handling ---#
  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster.", call. = FALSE)
  }

  if (!is.numeric(nStrata)) {
    stop("'nStrata' must be type numeric.", call. = FALSE)
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical.", call. = FALSE)
  }

  if (!is.numeric(samp)) {
    stop("'samp' must be type numeric.", call. = FALSE)
  }

  if (!is.logical(details)) {
    stop("'details' must be type logical.", call. = FALSE)
  }

  #--- if there is only 1 band in mraster use it as default ---#

  if (terra::nlyr(mraster) > 1) {
    stop("Multiple layers detected in 'mraster'. Please define a singular band to stratify.", call. = FALSE)
  }

  if (!is.null(mraster2)) {
    if (!inherits(mraster2, "SpatRaster")) {
      stop("'mraster2' must be type SpatRaster.", call. = FALSE)
    }

    if (isFALSE(terra::compareGeom(mraster, mraster2, stopOnError = FALSE))) {
      stop("Extents of 'mraster' and 'mraster2' do not match.", call. = FALSE)
    }
    
    if (isFALSE(terra::compareGeom(mraster, mraster2, stopOnError = FALSE, ext = FALSE, res = TRUE))) {
      stop("Spatial resolutions of 'mraster' and 'mraster2' do not match.", call. = FALSE)
    }

    #--- if there is only 1 band in mraster2 use it as default ---#

    if (terra::nlyr(mraster2) > 1) {
      stop("Multiple layers detected in 'mraster2'. Please define a singular band to stratify.", call. = FALSE)
    }

    #--- ensure names of rasters are not the same ---#

    if (identical(names(mraster), names(mraster2))) {
      stop("'mraster' and 'mraster2' must have different metric names.", call. = FALSE)
    }
    
    if (!is.numeric(nStrata2)) {
      stop("'nStrata2' must be type numeric.", call. = FALSE)
    }
  }

  if (is.null(mraster2)) {
    vals <- terra::values(mraster)

    vals[!is.finite(vals)] <- NA

    #--- Determine index of each cell so to map values correctly without NA ---#

    idx <- !is.na(vals)

    #--- Remove NA / NaN / Inf values ---#

    df <- vals %>%
      as.data.frame() %>%
      stats::na.omit()

    metric <- as.character(names(df))

    metric <- ggplot2::ensym(metric)

    #--- Split mraster distribution in to number specified by 'breaks' ---#

    dfc <- df %>%
      dplyr::mutate(class = dplyr::ntile(!!metric, nStrata))
    
    #--- generate lookup table for details output ---#
    lookUp <- dfc %>%
      dplyr::group_by(class) %>%
      dplyr::summarize(min_mraster = min(!!metric),
                    max_mraster = max(!!metric))

    #--- convert back to original mraster extent ---#

    suppressWarnings(vals[, 1][idx] <- dfc$class)

    #--- set newly stratified values ---#

    rout <- terra::setValues(mraster[[1]], vals[, 1])
    names(rout) <- "strata"

    #--- if mraster2 is provided ---#
  } else {

    #--- Extract values from mraster ---#

    vals <- c(mraster, mraster2) %>%
      terra::values()

    vals[!is.finite(vals)] <- NA

    #--- Determine index of each cell so to map values correctly without NA ---#

    idx <- is.finite(vals[, 1]) & is.finite(vals[, 2])

    #--- Remove NA / NaN / Inf values ---#

    df <- vals %>%
      as.data.frame() %>%
      stats::na.omit()

    nm <- as.character(names(df))

    metric <- nm[1]
    metric2 <- nm[2]

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
    
    #--- generate lookup table for details output ---#
    lookUp <- dfc %>%
      dplyr::group_by(class) %>%
      dplyr::summarize(min_mraster = min(!!metric),
                       max_mraster = max(!!metric),
                       min_mraster2 = min(!!metric2),
                       max_mraster2 = max(!!metric2))

    #--- convert back to original mraster extent ---#

    suppressWarnings(vals[, 1][idx] <- dfc$class)

    #--- set newly stratified values ---#

    rout <- terra::setValues(mraster[[1]], vals[, 1])
    names(rout) <- "strata"
  }

  if (isTRUE(plot)) {
    if (is.null(mraster2)) {

      #--- output histogram with quantile breaks ---#

      #--- determine numeric breaks ---#

      breaks <- dfc %>%
        dplyr::group_by(class) %>%
        dplyr::summarize(breaks = max(!!metric)) %>%
        dplyr::select(breaks) %>%
        as.data.frame()

      breaks <- breaks[1:(nStrata - 1), ]

      df.p <- dfc %>%
        dplyr::select(!!metric)

      #--- plot histogram of metric with associated break lines ---#

      p <- ggplot2::ggplot(df.p, ggplot2::aes(!!metric)) +
        ggplot2::geom_histogram() +
        ggplot2::geom_vline(xintercept = breaks, linetype = "dashed") +
        ggplot2::ggtitle(paste0(metric, " histogram with defined breaks"))

      terra::plot(rout, main = "Classes")
    } else {
      terra::plot(rout, main = "Classes")
    }
  }

  #--- write file to disc ---#

  if (!is.null(filename)) {
    terra::writeRaster(x = rout, filename = filename, overwrite = overwrite)
    message("Output raster written to disc.")
  }

  #--- Output based on 'details' to return raster alone or list with details ---#

  if (isTRUE(details)) {
    if (!is.null(mraster2)) {
      if (samp > 1 | samp < 0) {
        stop("'samp' must be > 0 <= 1.", call. = FALSE)
      }

      #--- create classplot summary ---#

      coordsgrps <- dfc %>%
        dplyr::group_by(class) %>%
        dplyr::arrange(class) %>%
        stats::na.omit() %>%
        tidyr::nest() %>%
        dplyr::ungroup()

      p <- classPlot(as.data.frame(dfc),
        coordsgrps = coordsgrps,
        mraster = !!metric,
        mraster2 = !!metric2,
        samp
      )
    }

    #--- output metrics details along with stratification raster ---#

    if (isTRUE(plot)){
    out <- list(
      details = lookUp,
      raster = rout,
      plot = p
    )
    } else {
      out <- list(
        details = lookUp,
        raster = rout
      )
    }

    return(out)
  } else {

    #--- just output raster ---#

    return(rout)
  }
}
