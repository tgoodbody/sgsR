#' Optimum sample breaks stratification
#'
#' @description Stratify metrics raster using optimum sample breaks algorithm
#' @family stratify functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams strat_breaks
#'
#' @param metric Character. Name of metric to be used for stratification
#' @param nStrata Numeric. Number of desired output strata.
#' @param nSamp Numeric. Number of desired samples - used within
#' OSB algorithm to help determine break points.
#' @param subset - Numeric. Value between 0 and 1 (default)
#' denoting proportion of data to use to determine break points
#'
#' @importFrom magrittr %>%
#' @importFrom methods is
#'
#' @return output stratification \code{spatRaster}, or a list when \code{details = TRUE}.
#'
#' @export

strat_osb <- function(mraster,
                      metric = NULL,
                      nStrata,
                      nSamp,
                      subset = 1,
                      plot = FALSE,
                      details = FALSE,
                      filename = NULL,
                      overwrite = FALSE,
                      ...) {

  #--- check for required packages ---#
  if (!requireNamespace("stratifyR", quietly = TRUE)) {
    stop("Package \"stratifyR\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  #--- Set global vars ---#

  from <- NULL

  #--- Error management ---#

  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  }

  if (!metric %in% names(mraster)) {
    stop(paste0("mraster does not have a variable named ", metric))
  }

  if (!is.numeric(nStrata)) {
    stop("'nStrata' must be type numeric")
  }

  if (!is.numeric(nSamp)) {
    stop("'nSamp' must be type numeric")
  }

  if (!is.numeric(subset)) {
    stop("'subset' must be type numeric")
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical")
  }

  if (!is.logical(details)) {
    stop("'details' must be type logical")
  }

  #--- if there is only 1 metric in the raster use it as default ---#

  if (terra::nlyr(mraster) == 1) {
    rastermetric <- mraster
  } else {
    if (is.null(metric)) {
      stop(" multiple layers detected in 'mraster'. Please define a 'metric' to stratify")
    }

    if (!is.character(metric)) {
      stop("'metric' must be type character")
    }

    #--- extract mraster metric ---#

    rastermetric <- terra::subset(mraster, metric)
  }

  #--- Perform OSB ---#
  #--- determine whether data should be subset prior to OSB calculation to save processing time ---#

  if (isTRUE(subset)) {
    if (subset > 1 | subset < 0) {
      stop("'subset' must be between 0 and 1")
    }

    message(paste0("'subset' was specified. Taking ", subset * 100, "% of available pixels to determine OSB"))

    #--- Extract values from mraster removing any NA/INF/NaN ---#

    OSB <- perform_osb_sample(rastermetric, nStrata, nSamp, subset)
  } else {
    if (terra::ncell(rastermetric) > 100000) {
      message("The raster you are using has over 100,000 cells. Consider using 'subset' to improve processing times.")
    }

    #--- Extract values from raster removing any NA/INF/NaN ---#

    OSB <- perform_osb(rastermetric, nStrata, nSamp)
  }

  #--- reclassify values based on breaks ---#

  breaks <- data.frame(from = c(-Inf, OSB[[2]]$OSB[1:(nStrata - 1)], Inf)) %>%
    dplyr::mutate(
      to = dplyr::lead(from),
      becomes = seq(1:length(from))
    ) %>%
    stats::na.omit() %>%
    as.matrix()

  rcl <- terra::classify(rastermetric, breaks)
  names(rcl) <- "strata"

  if (isTRUE(plot)) {
    data <- as.data.frame(OSB[[1]])
    names(data) <- "metric"

    #--- plot histogram of metric with associated break lines ---#

    p1 <- ggplot2::ggplot(data, ggplot2::aes(metric)) +
      ggplot2::geom_histogram() +
      ggplot2::geom_vline(xintercept = OSB[[2]]$OSB, linetype = "dashed") +
      ggplot2::ggtitle("Metric histogram with OSB break lines")

    print(p1)

    #--- set colour palette ---#

    ncols <- nStrata
    col <- RColorBrewer::brewer.pal(ncols, "Set3")

    terra::plot(rcl, main = "OSB breaks", col = col, type = "classes")
  }

  #--- write file to disc ---#

  if (!is.null(filename)) {
    terra::writeRaster(rcl, filename, overwrite = overwrite, ...)
  }

  #--- Output based on 'details' to return raster alone or list with details ---#

  if (isTRUE(details)) {

    #--- output OSB break points raster with associated breaks ---#

    breaks_rcl <- list(details = OSB[[2]]$OSB, raster = rcl)

    return(breaks_rcl)
  } else {

    #--- just output raster ---#

    return(rcl)
  }
}

perform_osb_sample <- function(rastermetric, nStrata, nSamp, subset) {
  vals <- rastermetric %>%
    terra::values(dataframe = TRUE) %>%
    dplyr::filter(stats::complete.cases(.)) %>%
    dplyr::slice_sample(prop = subset) %>%
    dplyr::pull()

  OSB_result <- vals %>%
    stratifyR::strata.data(h = nStrata, nSamp = nSamp)

  out <- list(vals, OSB_result)

  out
}

perform_osb <- function(rastermetric, nStrata, nSamp) {
  vals <- rastermetric %>%
    terra::values(dataframe = TRUE) %>%
    dplyr::filter(stats::complete.cases(.)) %>%
    dplyr::pull()

  OSB_result <- vals %>%
    stratifyR::strata.data(h = nStrata, nSamp = nSamp)

  out <- list(vals, OSB_result)

  out
}
