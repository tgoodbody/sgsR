#' Extract stratification raster strata to existing samples
#' @family extract functions
#'
#' @inheritParams  sample_srs
#' @inheritParams  extract_metrics
#' @param existing sf or data.frame.  Existing plot network.
#' 
#' @return An sf or data.frame object of samples with strata attributes
#' 
#' @export

extract_existing <- function(sraster,
                             existing,
                             data.frame = FALSE){
  
  #--- Error management ---#
  
  if (!inherits(sraster, "SpatRaster") )
    stop("'sraster' must be type SpatRaster", call. = FALSE)
  
  if (any(! c("strata") %in% names(sraster)) )
    stop("'sraster' must have a layer named 'strata'")
  
  if (!inherits(existing,"sf") && inherits(sf::st_geometry(existing),"sfc_POINT"))
    stop("'existing' must be an 'sf' object of type 'sfc_POINT' geometry")
  
  if (!is(existing, "data.frame"))
    stop("existing must be a data.frame")
  
  #--- if the existing plots are an sf object extract coordinates ---#
  
  if (is(existing, "sf")){
    
    #--- Convert to spatVector to enable extraction of strata values ---#
    
    existing <- sf::st_coordinates(existing)
    
  } else {
  
    if (any(! c("X", "Y") %in% colnames(existing)) ) {
    
      #--- if coordinate column names are lowercase change them to uppercase to match requirements ---#
      
      if (any(c("x", "y") %in% colnames(existing))) {
        
        existing <- existing %>%
          dplyr::rename(X = x,
                 Y = y)
        
        message("Column coordinate names are lowercase - converting to uppercase")
        
      } else {
      
        #--- if no x/y columns are present stop ---#  
        
        stop("'existing' must have columns named 'X' and 'Y'")
        
      }
      
    }
    
  }
  
  #--- extract values from the sraster dataset ---#
  
  strata_vals <- terra::extract(sraster,existing)
  
  #--- bind values and coordinates ---#
  
  existing_strata <- cbind(existing, strata_vals)
  
  #--- select only coordinate and strata values ---#
  
  existing_strata <- existing_strata %>% 
    dplyr::select(X,Y,strata)
  
  #--- output either data.frame or sf object ---#
  
  if (isTRUE(data.frame)) {
    
    #--- return data.frame ---#
    return(existing_strata)
    
  } else {
    
    #--- convert coordinates to a sf object ---#
    
    samples <- existing_strata %>%
      as.data.frame() %>%
      sf::st_as_sf(., coords = c("X", "Y"))
    
    #--- assign sraster crs to spatial points object ---#
    
    sf::st_crs(samples) <- terra::crs(sraster)
    
    #--- return sf object ---#
    return(samples)
    
  }
  
}