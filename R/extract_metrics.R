#' Extract metric raster attributes to samples
#' @family extract functions
#'
#' @inheritParams sample_srs
#' @inheritParams strat_kmeans
#' @param samples sf. Samples resulting from sample_* functions.
#' @param data.frame Logical. If true outputs as data.frame
#' 
#' @return An sf or data.frame object of samples with associated raster cell attributes
#' 
#' @export

extract_metrics <- function(mraster,
                            samples,
                            data.frame = FALSE){

  #--- Error management ---#
  
  if (!inherits(mraster, "SpatRaster"))
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  
  if (!inherits(samples, "sf"))
    stop("'samples' must be an 'sf' object", call. = FALSE)
  
  #--- Extract coordinates from samples ---#
  
  xy <- sf::st_coordinates(samples)
  
  vals <- terra::extract(mraster,xy)
  
  #--- extract other attributes from sampling and remove geometry attribute ---#
  
  samp_mets <- as.data.frame(samples)
  
  samp_mets <-  dplyr::select(samp_mets, -geometry)
  
  #--- bind values and coordinates ---#
  samples <- cbind(xy, samp_mets, vals)
  
  if (isTRUE(data.frame)) {
    
    #--- return data.frame ---#
    return(samples)
    
  } else {
    
    #--- convert coordinates to a sf object ---#
    
    samples <- samples %>%
      as.data.frame() %>%
      sf::st_as_sf(., coords = c("X", "Y"))
    
    #--- assign mraster crs to spatial points object ---#
    
    sf::st_crs(samples) <- terra::crs(mraster)
    
    #--- return sf object ---#
    return(samples)
    
  }

}