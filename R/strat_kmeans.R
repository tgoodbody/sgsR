## need to filter initial als metrics to a select set to be used for structurally guided sampling
##

strat_kmeans <- function(raster,
                         k,
                         iter.max = 500,
                         algorithm = "Lloyd",
                         center = TRUE,
                         scale = TRUE,
                         plot = TRUE
)
{

  #--- Error management ---#
  if (!inherits(raster,"SpatRaster"))
    stop("'raster_in' must be type SpatRaster", call. = FALSE)

  if (!is.numeric(k))
    stop("'k' must be type numeric")

  if (!is.numeric(iter.max))
    stop("'iter.max' must be type numeric")

  if (!is.character(algorithm))
    stop("'algorithm' must be type character")

  if (algorithm != 'Lloyd' && algorithm != 'MacQueen')
    stop("Unknown algorithm '", algorithm,"' selected. Please use 'Lloyd' (default) or 'MacQueen'")

  if (!is.logical(center))
    stop("'center' must be type logical")

  if (!is.logical(scale))
    stop("'scale' must be type logical")

  if (!is.logical(plot))
    stop("'plot' must be type logical")


  #--- Extract values from raster ---#
  vals <- terra::values(raster)

  #--- Determine index of each cell so to map values correctly without NA ---#
  vals[!is.finite(vals)] <- NA

  idx <- !is.na(vals)

  #--- conduct unsupervised k-means with center/scale parameters based on algorithm ---#
  if (algorithm == "Lloyd")
    {
    if (center == TRUE & scale == TRUE)
      {
      message("K-means being performed on ",nlyr(raster)," layers with ",k," centers. Data have been centered and scaled")

      km_clust <- stats::kmeans(scale(vals[idx],center = TRUE, scale=TRUE), centers=k, iter.max = iter.max, algorithm = algorithm)
    } else if (center == TRUE & scale == FALSE) {
      message("K-means being performed on ",nlyr(raster)," layers with ",k," centers. Data have been centered but NOT scaled.")

      km_clust <- stats::kmeans(scale(vals[idx],center = TRUE, scale=FALSE), centers=k, iter.max = iter.max, algorithm = algorithm)
    } else if (center == FALSE & scale == TRUE) {
      message("K-means being performed on ",nlyr(raster)," layers with ",k," centers. Data have been scaled but NOT centered.")

      km_clust <- stats::kmeans(scale(vals[idx],center = FALSE, scale=TRUE), centers=k, iter.max = iter.max, algorithm = algorithm)
    } else if (center == FALSE & scale == FALSE) {
      message("K-means being performed on ",nlyr(raster)," layers with ",k," centers, Data have not been centered or scaled.")

      km_clust <- stats::kmeans(scale(vals[idx],center = FALSE, scale=FALSE), centers=k, iter.max = iter.max, algorithm = algorithm)
    }
  }
  #--- conduct unsupervised k-means with center/scale parameters based on algorithm ---#
  else if (algorithm == "MacQueen")
    {
    if (center == TRUE & scale == TRUE){
      message("K-means being performed on ",nlyr(raster)," layers with ",k," centers. Data have been centered and scaled")

      km_clust <- stats::kmeans(scale(vals[idx],center = TRUE, scale=TRUE), centers=k, iter.max = iter.max, algorithm = algorithm)
    } else if (center == TRUE & scale == FALSE) {
      message("K-means being performed on ",nlyr(raster)," layers with ",k," centers. Data have been centered but NOT scaled.")

      km_clust <- stats::kmeans(scale(vals[idx],center = TRUE, scale=FALSE), centers=k, iter.max = iter.max, algorithm = algorithm)
    } else if (center == FALSE & scale == TRUE) {
      message("K-means being performed on ",nlyr(raster)," layers with ",k," centers. Data have been scaled but NOT centered.")

      km_clust <- stats::kmeans(scale(vals[idx],center = FALSE, scale=TRUE), centers=k, iter.max = iter.max, algorithm = algorithm)
    } else if (center == FALSE & scale == FALSE) {
      message("K-means being performed on ",nlyr(raster)," layers with ",k," centers, Data have not been centered or scaled.")

      km_clust <- stats::kmeans(scale(vals[idx],center = FALSE, scale=FALSE), centers=k, iter.max = iter.max, algorithm = algorithm)
    }
  }

  #--- convert k-means values back to original raster extent ---#
  vals[idx] <- km_clust$cluster

  kmv <- terra::setValues(raster[[1]],vals)
  names(kmv) <- "class"

  #--- assign output raster to k-means info ---#
  km_clust$raster <- kmv

    #--- plot if requested ---#
  if (plot == TRUE){

    #--- make plot using diverging colour palette ---#
    ncols <- k
    col = brewer.pal(ncols, "Set3")

    terra::plot(kmv, main = 'K-means clusters', col=col,type="classes")

    #--- return k-means object with combined raster ---#
    return(km_clust)

  }
    #--- return k-means object with combined raster ---#
    return(km_clust)

}
