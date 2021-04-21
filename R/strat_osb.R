#' Stratify raster using optimum sample breaks algorithm
#' @family stratify functions
#'
#' @inheritParams strat_kmeans
#' @param raster spatRaster. Multiband ALS metrics raster.
#' @param metric Character. Name of metric to be used for stratification
#' @param h Numeric. Number of desired output strata.
#' @param n Numeric. Number of desired samples - used within OSB algorithm to help determine break points.
#' @param subset - Logical. Determines whether a subset of data should be used to help determine break points
#'
#' @return list where \code{breaks} are the breaks defined by the OSB algorithm and \code{raster} is the output stratification spatRaster
#' 
#' @export

strat_osb <- function(raster,
                      metric,
                      h,
                      n,
                      subset = TRUE,
                      plot = FALSE)

{
  #--- Error management ---#
  if (!inherits(raster,"SpatRaster"))
    stop("'raster' must be type SpatRaster", call. = FALSE)
  
  if (!is.character(metric))
    stop("'metric' must be type character")
  
  if(! metric %in% names(raster))
    stop(paste0("raster does not have a variable named ",metric))

  if (!is.numeric(h))
    stop("'h' must be type numeric")

  if (!is.numeric(n))
    stop("'n' must be type numeric")

  if (!is.logical(subset))
    stop("'subset' must be type logical")

  if (!is.logical(plot))
    stop("'plot' must be type logical")
  
  #--- extract raster metric ---#
  rastermetric <- terra::subset(raster,metric)

  #--- Perform OSB ---#
  #--- determine whether data should be subset prior to OSB calculation to save processing time ---#
  if (subset == TRUE) {
    message("'subset = TRUE' : taking 50% of available pixels to determine OSB")

    #--- Extract values from raster removing any NA/INF/NaN ---#
    OSB <- perform_osb_sample(rastermetric, h, n)

  } else {

    if (ncell(rastermetric) > 100000)
      message("The raster you are using has over 100,000 cells. Consider using 'subset = TRUE' to improve processing times.")

    #--- Extract values from raster removing any NA/INF/NaN ---#
    OSB <- perform_osb(rastermetric, h, n)

  }

  #--- reclassify values based on breaks ---#

  breaks <- data.frame(from=c(-Inf,OSB[[2]]$OSB[1:(h - 1)],Inf)) %>%
    mutate(to = lead(from),
           becomes = seq(1:length(from))) %>%
    na.omit() %>%
    as.matrix()

  rcl <- terra::classify(rastermetric,breaks)
  names(rcl) <- "strata"

  if (plot == TRUE){

    data <- as.data.frame(OSB[[1]])
    names(data) <- "metric"

    #--- plot histogram of metric with associated break lines ---#

    p1 <- ggplot(data,aes(metric)) +
      geom_histogram() +
      geom_vline(xintercept = OSB[[2]]$OSB, linetype = "dashed") +
      ggtitle("Metric histogram with OSB break lines")

    print(p1)

    #--- set colour palette ---#
    ncols <- h
    col = brewer.pal(ncols, "Set3")

    terra::plot(rcl, main = 'OSB breaks', col=col,type="classes")

  }

  #--- output OSB break points raster with associated breaks ---#
  breaks_rcl <- list(breaks = OSB[[2]]$OSB, raster = rcl)

  return(breaks_rcl)

}

perform_osb_sample <- function(rastermetric, h, n){
  vals <- rastermetric %>%
    terra::values(dataframe=TRUE) %>%
    filter(complete.cases(.)) %>%
    slice_sample(prop = 0.5) %>%
    pull()

  OSB_result <- vals %>%
    strata.data(h = (h), n = n)

  out <- list(vals,OSB_result)

  out
}

perform_osb <- function(rastermetric, h, n){
  vals <- rastermetric %>%
    terra::values(dataframe=TRUE) %>%
    filter(complete.cases(.)) %>%
    pull()

  OSB_result <- vals %>%
    strata.data(h = (h), n = n)

  out <- list(vals,OSB_result)

  out
}

