

strat_osb <- function(band,
                      h,
                      n,
                      subset = TRUE,
                      plot = TRUE,
                      ...
)

{
  #--- Error management ---#
  if (!inherits(band,"SpatRaster"))
    stop("all specified bands must be type SpatRaster", call. = FALSE)

  if (!is.numeric(h))
    stop("'h' must be type numeric")

  if (!is.numeric(n))
    stop("'n' must be type numeric")

  if (!is.logical(subset))
    stop("'subset' must be type logical")

  if (!is.logical(plot))
    stop("'plot' must be type logical")

  #--- Perform OSB ---#
  #--- determine whether data should be subset prior to OSB calculation to save processing time ---#
  if (subset == TRUE) {
    message("'subset = TRUE' : taking 50% of available pixels to determine OSB")

    #--- Extract values from raster removing any NA/INF/NaN ---#
    OSB <- perform_osb_sample(band, h, n)

  } else {

    if (ncell(band) > 100000)
      message("The raster you are using has over 100,000 cells. Consider using 'subset = TRUE' to improve processing times.")

    #--- Extract values from raster removing any NA/INF/NaN ---#
    OSB <- perform_osb(band, h, n)

  }

  #--- reclassify values based on breaks ---#

  breaks <- data.frame(from=c(-Inf,OSB[[2]]$OSB)) %>%
    mutate(to = lead(from),
           becomes = seq(1:length(from))) %>%
    na.omit() %>%
    as.matrix()

  rcl <- terra::classify(band,breaks,include.lowest=TRUE)
  names(rcl) <- "class"

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
  breaks_rcl <- list(OSB[[2]]$OSB,rcl)

  return(breaks_rcl)

}

perform_osb_sample <- function(band,h,n){
  vals <- band %>%
    terra::values(dataframe=TRUE) %>%
    filter(complete.cases(.)) %>%
    slice_sample(prop = 0.5) %>%
    pull()

  OSB_result <- vals %>%
    strata.data(h = (h), n = n)

  out <- list(vals,OSB_result)

  out
}

perform_osb <- function(band,h,n){
  vals <- band %>%
    terra::values(dataframe=TRUE) %>%
    filter(complete.cases(.)) %>%
    pull()

  OSB_result <- vals %>%
    strata.data(h = (h), n = n)

  out <- list(vals,OSB_result)

  out
}

