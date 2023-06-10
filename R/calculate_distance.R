#' Distance to access layer
#'
#' @description Per pixel distance to nearest access vector. Intended to be used as a `cost` constraint
#' within the \code{\link{sample_clhs}} function
#'
#' @family calculate functions
#'
#' @inheritParams sample_srs
#' @inheritParams strat_breaks
#'
#' @param raster spatRaster. Raster to be used to calculate pixel level distance to access layer.
#' @param slope Logical. Calculate slope distance instead of geographic distance. \code{raster} needs 
#' to be a digital terrain model for this to make sense.
#'
#' @return Input raster with \code{dist2access} layer appended.
#'
#' @examples
#' \dontrun{
#' #--- Load raster and access files ---#
#' r <- system.file("extdata", "mraster_small.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' a <- system.file("extdata", "access.shp", package = "sgsR")
#' ac <- sf::st_read(a)
#'
#' calculate_distance(
#'   raster = mr,
#'   access = ac,
#' )
#' }
#'
#' @author Tristan R.H. Goodbody
#'
#' @export

calculate_distance <- function(raster,
                               access,
                               slope = FALSE,
                               plot = FALSE,
                               filename = NULL,
                               overwrite = FALSE) {
  #--- error handling ---#

  if (!inherits(raster, "SpatRaster")) {
    stop("'raster' must be type SpatRaster.", call. = FALSE)
  }

  if (!inherits(access, "sf")) {
    stop("'access' must be an 'sf' object.", call. = FALSE)
  }
  
  if (!is.logical(slope)) {
    stop("'slope' must be type logical.", call. = FALSE)
  }

  if (!inherits(sf::st_geometry(access), "sfc_MULTILINESTRING") && !inherits(sf::st_geometry(access), "sfc_LINESTRING")) {
    stop("'access' geometry type must be 'LINESTRING' or 'MULTILINESTRING'.", call. = FALSE)
  }

  #--- load access ---#

  access <- terra::vect(access)

  if (terra::crs(access) != terra::crs(raster)) terra::crs(access) <- terra::crs(raster)

  #--- use first layer from raster and access to determine distance from each pixel ---#

  message("calculating per pixel distance to provided access layer")
  
  if(isFALSE(slope)){
    
    dist2access <- terra::distance(raster[[1]], access)
    
  } else {
    
    #--- from @spono ---#
    
    slopeDist <- terra::terrain(raster[[1]], v = "slope", unit = "radians", neighbors = 8)
    
    resolution <- terra::res(slopeDist)[1]
    
    slopeDist <-  resolution / cos(slopeDist)
    
    access$value <- 0
    
    access_rast <- terra::rasterize(access, slopeDist, field = "value", background = 1) # values=0
    
    slopeDist <- slopeDist * access_rast
    
    dist2access <- terra::costDist(slopeDist, target = 0, scale = resolution, maxiter = 50)
    
  }

  #--- append dist2access layer to raster ---#

  raster$dist2access <- dist2access

  if (isTRUE(plot)) {
    terra::plot(dist2access)
    suppressWarnings(terra::plot(access, add = T))
  }

  write_raster(raster = raster, filename = filename, overwrite = overwrite)

  return(raster)
}
