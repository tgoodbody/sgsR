# raster = spatRaster. Multiband ALS metrics raster.
# ncp = Character. Number of components to create.
# b1 = Numeric. Number of desired strata for first principal component.
# b2 = Numeric. Number of desired strata for second principal component.
# scale - Logical. Determines whether centering and scaling of data should be conducted prior to principal component analysis.
# plot = Logical. Plots output strata raster and visualized strata with boundary dividers.
# samp = Numeric. Determines proportion of cells to plot for strata visualization. Lower values reduce processing time.

## output is a list where '$breaks' are the breaks defined by the OSB algorithm and '$raster' is the output stratification spatRaster


strat_osb <- function(raster,
                      metric,
                      h,
                      n,
                      subset = TRUE,
                      plot = TRUE,
                      ...
)

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

  breaks <- data.frame(from=c(-Inf,OSB[[2]]$OSB)) %>%
    mutate(to = lead(from),
           becomes = seq(1:length(from))) %>%
    na.omit() %>%
    as.matrix()

  rcl <- terra::classify(rastermetric,breaks,include.lowest=TRUE)
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

perform_osb_sample <- function(rastermetric,h,n){
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

perform_osb <- function(rastermetric,h,n){
  vals <- rastermetric %>%
    terra::values(dataframe=TRUE) %>%
    filter(complete.cases(.)) %>%
    pull()

  OSB_result <- vals %>%
    strata.data(h = (h), n = n)

  out <- list(vals,OSB_result)

  out
}

