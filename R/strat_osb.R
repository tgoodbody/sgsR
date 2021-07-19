#' Determine optimum sample boundaries
#'
#' @description Determine optimum sample boundaries algorithm of univariate populations
#' using the \code{\link[stratifyR]{strata.data}} algorithm.
#'
#' @family stratify functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams strat_breaks
#'
#' @param metric Character. Name of metric to be used for stratification.
#' @param nStrata Numeric. Number of desired output strata.
#' @param nSamp Numeric. Number of desired samples - used within
#' OSB algorithm to help determine break points.
#' @param subset Numeric. Value between 0 and 1 (default)
#' denoting proportion of data to use to determine break points.
#'
#'
#' @references
#' Khan, E. A., Khan, M. G. M., & Ahsan, M. J. (2002). Optimum Stratification:
#' A Mathematical Programming Approach. Calcutta Statistical Association Bulletin,
#' 52(1–4), 323–334. https://doi.org/10.1177/0008068320020518
#'
#' Khan, M. G. M., Nand, N., & Ahmad, N. (2008). Determining the optimum strata
#' boundary points using dynamic programming. Survey methodology, 34(2), 205-214.
#'
#' M.G.M. Khan, K.G. Reddy & D.K. Rao (2015) Designing stratified sampling in economic
#' and business surveys, Journal of Applied Statistics, 42:10, 2080-2099,
#' DOI: 10.1080/02664763.2015.1018674
#'
#' @return Returns an output stratification \code{spatRaster} or a list when \code{details = TRUE}.
#'
#' When a list is returned:
#' \enumerate{
#' \item \code{details} is a list output of the \code{\link[stratifyR]{strata.data}} function where
#' \code{OSB} are the optimum stratum boundaries and \code{nh} are the optimum sample sizes
#' for each stratum.
#' \item \code{osb} vector of optimum stratum boundaries.
#' \item \code{breaks} matrix associated metric and strata break values.
#' \item \code{raster} is a stratified \code{spatRaster} based on \code{OSB}.
#' }
#'
#'
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
      message("Consider using 'subset' to improve processing times.")
    }
    
    #--- Extract values from raster removing any NA/INF/NaN ---#
    
    OSB <- perform_osb(rastermetric, nStrata, nSamp)
  }
  
  #--- reclassify values based on breaks ---#
  
  breaks <- data.frame(from = c(-Inf, OSB[[2]]$OSB[1:(nStrata - 1)], Inf)) %>%
    dplyr::mutate(
      to = dplyr::lead(from),
      becomes = seq(1:length(from))) %>%
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
    
    breaks_rcl <- list(details = OSB, osb = OSB[[2]]$OSB, breaks = breaks, raster = rcl)
    
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
    stratifyR::strata.data(h = nStrata, n = nSamp)
  
  out <- list(vals, OSB_result)
  
  out
}

perform_osb <- function(rastermetric, nStrata, nSamp) {
  vals <- rastermetric %>%
    terra::values(dataframe = TRUE) %>%
    dplyr::filter(stats::complete.cases(.)) %>%
    dplyr::pull()
  
  OSB_result <- vals %>%
    stratifyR::strata.data(h = nStrata, n = nSamp)
  
  out <- list(vals, OSB_result)
  
  out
}