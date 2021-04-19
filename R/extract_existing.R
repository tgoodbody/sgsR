# raster = spatRaster. Stratification raster derived from sample_* functions.
# existing = sf or data.frame.  Existing plot network to be included within sample_* functions
# data.frame = Logical. If TRUE, the user receives a data.frame output instead of an sf object

extract_existing <- function(raster,
                            existing,
                            data.frame = FALSE){
  
  #--- Error management ---#
  
  if (!inherits(raster, "SpatRaster") )
    stop("'raster' must be type SpatRaster", call. = FALSE)
  
  if (any(! c("strata") %in% names(raster)) )
    stop("'raster' must have a layer named 'strata'")
  
  if (!inherits(existing,"sf") && inherits(sf::st_geometry(existing),"sfc_POINT"))
    stop("'sf' object must be of type 'sfc_POINT' geometry")
  
  if (!is(existing, "data.frame"))
    stop("existing must be a data.frame")
  
  #--- if the existing plots are an sf object extract coordinates ---#
  
  if (is(existing, "sf")){
    
    #--- Convert to spatVector to enable extraction of strata values ---#
    
    existing <- st_coordinates(existing)
    
  } else {
  
    if (any(! c("X", "Y") %in% colnames(existing)) ) {
    
      #--- if coordinate column names are lowercase change them to uppercase to match requirements ---#
      
      if (any(c("x", "y") %in% colnames(existing))) {
        
        existing <- existing %>%
          rename(X = x,
                 Y = y)
        
        message("Column coordinate names are lowercase - converting to uppercase")
        
      } else {
      
        #--- if no x/y columns are present stop ---#  
        
        stop("'existing' must have columns named 'X' and 'Y'")
        
      }
      
    }
    
  }
  
  #--- extract values from the raster dataset ---#
  
  strata_vals <- terra::extract(raster,existing)
  
  #--- bind values and coordinates ---#
  
  existing_strata <- cbind(existing, strata_vals)
  
  #--- select only coordinate and strata values ---#
  
  existing_strata <- existing_strata %>% 
    dplyr::select(X,Y,strata)
  
  #--- output either data.frame or sf object ---#
  
  if (data.frame == TRUE) {
    
    #--- return data.frame ---#
    return(existing_strata)
    
  } else {
    
    #--- convert coordinates to a sf object ---#
    
    samples <- existing_strata %>%
      as.data.frame() %>%
      st_as_sf(., coords = c("X", "Y"))
    
    #--- assign raster crs to spatial points object ---#
    
    st_crs(samples) <- crs(raster)
    
    #--- return sf object ---#
    return(samples)
    
  }
  
}