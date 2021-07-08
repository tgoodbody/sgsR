#' Optimum sample breaks stratification
#' 
#' @description Stratify metrics raster using optimum sample breaks algorithm
#' @family stratify functions
#'
#' @inheritParams strat_kmeans
#' @param metric Character. Name of metric to be used for stratification
#' @param nstrata Numeric. Number of desired output strata.
#' @param n Numeric. Number of desired samples - used within 
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
                      nstrata,
                      n,
                      subset = 1,
                      plot = FALSE,
                      details = FALSE)

{
  
  #--- Set global vars ---#
  
  from <- NULL
  
  #--- Error management ---#
  
  if (!inherits(mraster,"SpatRaster"))
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  
  if(! metric %in% names(mraster))
    stop(paste0("mraster does not have a variable named ",metric))

  if (!is.numeric(nstrata))
    stop("'nstrata' must be type numeric")

  if (!is.numeric(n))
    stop("'n' must be type numeric")

  if (!is.numeric(subset))
    stop("'subset' must be type numeric")

  if (!is.logical(plot))
    stop("'plot' must be type logical")
  
  if (!is.logical(details))
    stop("'details' must be type logical")
  
  #--- if there is only 1 metric in the raster use it as default ---#
  
  if (terra::nlyr(mraster) == 1){
    
    rastermetric <- mraster  
    
  } else {
    
    if (is.null(metric))
      stop(" multiple layers detected in 'mraster'. Please define a 'metric' to stratify")
    
    if (!is.character(metric))
      stop("'metric' must be type character")
  
    #--- extract mraster metric ---#
    
    rastermetric <- terra::subset(mraster,metric)
  
  }

  #--- Perform OSB ---#
  #--- determine whether data should be subset prior to OSB calculation to save processing time ---#
  
  if (isTRUE(subset)) {
    if (subset > 1 | subset < 0)
      stop("'subset' must be between 0 and 1")
    
    message(paste0("'subset' was specified. Taking ", subset * 100, "% of available pixels to determine OSB"))

    #--- Extract values from mraster removing any NA/INF/NaN ---#
    
    OSB <- perform_osb_sample(rastermetric, nstrata, n, subset)

  } else {

    if (terra::ncell(rastermetric) > 100000)
      message("The raster you are using has over 100,000 cells. Consider using 'subset' to improve processing times.")

    #--- Extract values from raster removing any NA/INF/NaN ---#
    
    OSB <- perform_osb(rastermetric, nstrata, n)

  }

  #--- reclassify values based on breaks ---#

  breaks <- data.frame(from = c(-Inf,OSB[[2]]$OSB[1:(nstrata - 1)],Inf)) %>%
    dplyr::mutate(to = dplyr::lead(from),
           becomes = seq(1:length(from))) %>%
    stats::na.omit() %>%
    as.matrix()

  rcl <- terra::classify(rastermetric,breaks)
  names(rcl) <- "strata"

  if (isTRUE(plot)){

    data <- as.data.frame(OSB[[1]])
    names(data) <- "metric"

    #--- plot histogram of metric with associated break lines ---#

    p1 <- ggplot2::ggplot(data, ggplot2::aes(metric)) +
      ggplot2::geom_histogram() +
      ggplot2::geom_vline(xintercept = OSB[[2]]$OSB, linetype = "dashed") +
      ggplot2::ggtitle("Metric histogram with OSB break lines")

    print(p1)

    #--- set colour palette ---#
    
    ncols <- nstrata
    col = RColorBrewer::brewer.pal(ncols, "Set3")

    terra::plot(rcl, main = 'OSB breaks', col=col,type="classes")

  }
  
  #--- Output based on 'details' to return raster alone or list with details ---#
  
  if ( isTRUE(details) ){
    
    #--- output OSB break points raster with associated breaks ---#
    
    breaks_rcl <- list(details = OSB[[2]]$OSB, raster = rcl)
    
    return(breaks_rcl)
    
    
  } else {
    
    #--- just output raster ---#
    
    return(rcl)
    
  }

}

perform_osb_sample <- function(rastermetric, nstrata, n, subset){
  vals <- rastermetric %>%
    terra::values(dataframe=TRUE) %>%
    dplyr::filter(stats::complete.cases(.)) %>%
    dplyr::slice_sample(prop = subset) %>%
    dplyr::pull()

  OSB_result <- vals %>%
    stratifyR::strata.data(h = nstrata, n = n)

  out <- list(vals,OSB_result)

  out
}

perform_osb <- function(rastermetric, nstrata, n){
  
  vals <- rastermetric %>%
    terra::values(dataframe=TRUE) %>%
    dplyr::filter(stats::complete.cases(.)) %>%
    dplyr::pull()

  OSB_result <- vals %>%
    stratifyR::strata.data(h = nstrata, n = n)

  out <- list(vals,OSB_result)

  out
}

