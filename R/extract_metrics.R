# raster = spatRaster. ALS metrics raster used to create the stratification rasters.
# samples = sf.  Samples resulting from sample_* functions.

extract_metrics <- function(raster,
                            samples,
                            data.frame = FALSE){

  #--- Error management ---#
  
  if (!inherits(raster, "SpatRaster"))
    stop("'raster' must be type SpatRaster", call. = FALSE)
  
  if (!inherits(samples, "sf"))
    stop("'samples' must be an 'sf' object", call. = FALSE)
  
  #--- Extract coordinates from samples ---#
  
  xy <- st_coordinates(samples)
  
  vals <- terra::extract(raster,xy)
  
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
      st_as_sf(., coords = c("X", "Y"))
    
    #--- assign raster crs to spatial points object ---#
    
    st_crs(samples) <- crs(raster)
    
    #--- return sf object ---#
    return(samples)
    
  }

}