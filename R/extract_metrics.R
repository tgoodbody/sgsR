#' Extract metric raster attributes to samples
#' @family extract functions
#'
#' @inheritParams  sample_srs
#' @param samples sf. Samples resulting from sample_* functions.
#' @param data.frame Logical. If true outputs as data.frame
#' 
#' @return An sf or data.frame object of samples with associated raster cell attributes
#' 
#' @export

extract_metrics <- function(sraster,
                            samples,
                            data.frame = FALSE){

  #--- Error management ---#
  
  if (!inherits(sraster, "SpatRaster"))
    stop("'sraster' must be type SpatRaster", call. = FALSE)
  
  if (!inherits(samples, "sf"))
    stop("'samples' must be an 'sf' object", call. = FALSE)
  
  #--- Extract coordinates from samples ---#
  
  xy <- st_coordinates(samples)
  
  vals <- terra::extract(sraster,xy)
  
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
    
    #--- assign sraster crs to spatial points object ---#
    
    st_crs(samples) <- crs(sraster)
    
    #--- return sf object ---#
    return(samples)
    
  }

}