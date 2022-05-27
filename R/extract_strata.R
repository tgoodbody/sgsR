#' Extract strata
#' 
#' @description Extract stratum values to existing samples
#'
#' @family extract functions
#' 
#' @inheritParams sample_systematic
#' 
#' @param sraster spatRaster. Stratification raster.
#' @param existing sf.  Existing plot network.
#' @param data.frame Logical. Output as data.frame if \code{TRUE}
#' 
#' @return An sf or data.frame object of samples with strata attribute
#' 
#' @examples 
#' #--- Load sraster ---#
#' r <- system.file("extdata", "kmeans.tif", package = "sgsR")
#' sr <- terra::rast(r)
#' 
#' #--- load existing samples ---#
#' e <- system.file("extdata", "existing.shp", package = "sgsR")
#' e <- sf::st_read(e)
#' 
#' extract_strata(sraster = sr,
#'                existing = e)
#' 
#' @author Tristan R.H. Goodbody
#' 
#' @export

extract_strata <- function(sraster,
                           existing,
                           data.frame = FALSE,
                           filename = NULL,
                           overwrite = FALSE) {
  
  #--- Set global vars ---#
  
  x <- y <- X <- Y <- strata <- NULL
  
  #--- Error management ---#
  
  if (!inherits(sraster, "SpatRaster")) {
    stop("'sraster' must be type SpatRaster", call. = FALSE)
  }
  
  if (any(!c("strata") %in% names(sraster))) {
    stop("'sraster' must have a layer named 'strata'", call. = FALSE)
  }
  
  if (!inherits(existing, "sf") && inherits(sf::st_geometry(existing), "sfc_POINT")) {
    stop("'existing' must be an 'sf' object of type 'sfc_POINT' geometry", call. = FALSE)
  }
  
  if (!is(existing, "data.frame")) {
    stop("existing must be a data.frame", call. = FALSE)
  }
  
  #--- if the existing plots are an sf object extract coordinates ---#
  
  if (is(existing, "sf")) {
    
    #--- Convert to spatVector to enable extraction of strata values ---#
    
    existing <- sf::st_coordinates(existing)
  } else {
    if (any(!c("X", "Y") %in% colnames(existing))) {
      
      #--- if coordinate column names are lowercase change them to uppercase to match requirements ---#
      
      if (any(c("x", "y") %in% colnames(existing))) {
        existing <- existing %>%
          dplyr::rename(
            X = x,
            Y = y
          )
        
        message("Column coordinate names are lowercase - converting to uppercase")
      } else {
        
        #--- if no x/y columns are present stop ---#
        
        stop("'existing' must have columns named 'X' and 'Y'")
      }
    }
  }
  
  #--- extract values from the sraster dataset ---#
  
  strata_vals <- terra::extract(sraster, existing)
  
  #--- bind values and coordinates ---#
  
  existing_strata <- cbind(existing, strata_vals)
  
  #--- select only coordinate and strata values ---#
  
  existing_strata <- existing_strata %>%
    dplyr::select(X, Y, strata)
  
  #--- if existing samples are not linked with a stratum ---#
  if(any(!complete.cases(existing_strata$strata))){
    
    nNA <- existing_strata %>%
      filter(!complete.cases(strata)) %>%
      tally() %>%
      pull()
    
    message(paste0(nNA," samples are located where strata values are NA."))
  }
  
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
    
    if (!is.null(filename)) {
      if (!is.logical(overwrite)) {
        stop("'overwrite' must be either TRUE or FALSE")
      }
      
      if (file.exists(filename) & isFALSE(overwrite)) {
        stop(paste0("'",filename, "' already exists and overwrite = FALSE"))
      }
      
      sf::st_write(samples, filename, delete_layer = overwrite)
    }
    
    #--- return sf object ---#
    return(samples)
  }
}
