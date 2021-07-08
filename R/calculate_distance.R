#' distance to access 
#' 
#' @description Per pixel distance to nearest access vector. Intended to be used as a `cost` constraint
#' within the \code{sample_clhs} function
#' 
#' @family calculate functions
#' 
#' @inheritParams sample_srs
#' 
#' @param raster spatRaster. Raster to be used to calculate pixel level distance to access layer.
#'
#' @importFrom magrittr %>%
#' @importFrom methods is
#' 
#' @return Input raster with `dist2access` layer appended.
#' 
#' @export

calculate_distance <- function(raster = NULL,
                               access = NULL,
                               plot = FALSE)
{
  
  #--- error handling ---#
  
  if (!inherits(raster, "SpatRaster"))
    stop("'raster' must be type SpatRaster", call. = FALSE)
  
  if (is.na(terra::crs(raster)))
    stop("'raster' does not have a coordinate system")
  
  if (!inherits(access,"sf"))
    stop("'access' must be an 'sf' object")
  
  if (!inherits(sf::st_geometry(access),"sfc_MULTILINESTRING"))
    stop("'access' geometry type must be 'sfc_MULTILINESTRING'")

  #--- load access ---#
  
  access <- terra::vect(access)
  
  #--- use first layer from raster and access to determine distance from each pixel ---#
  
  message("calculating per pixel distance to provided access layer")
  
  dist2access <- terra::distance(raster[[1]], access)
  
  #--- append dist2access layer to raster ---#
  
  raster$dist2access <- dist2access
  
  if(isTRUE(plot)){
    
    terra::plot(dist2access)
    suppressWarnings(terra::plot(access, add = T))
                     
  }
  
  return(raster)
  
}

