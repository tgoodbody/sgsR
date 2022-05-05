#' Extract metrics
#' 
#' @description Extract metric values to existing samples
#'
#' @family extract functions
#' 
#' @inheritParams sample_systematic
#' 
#' @param mraster spatRaster. Metrics Raster.
#' @param existing sf.  Existing plot network.
#' @param data.frame Logical. Output as data.frame if \code{TRUE}
#' 
#' @return An sf or data.frame object of samples with metrics attributes
#' 
#' @examples 
#' #--- Load mraster ---#
#' r <- system.file("extdata", "wall_metrics.tif", package = "sgsR")
#' mr <- terra::rast(r)
#' 
#' #' #--- load existing samples ---#
#' e <- system.file("extdata", "existing.shp", package = "sgsR")
#' e <- sf::st_read(e)
#' 
#' extract_metrics(mraster =  mr,
#'                 existing = e)
#' 
#' @author Tristan R.H. Goodbody
#' 
#' @export

extract_metrics <- function(mraster,
                            existing,
                            data.frame = FALSE,
                            filename = NULL,
                            overwrite = FALSE) {
  
  #--- Set global vars ---#
  
  x <- y <- X <- Y <- strata <- geometry <- NULL
  
  #--- Error management ---#
  
  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  }
  
  if (!inherits(existing, "sf")) {
    stop("'existing' must be an 'sf' object", call. = FALSE)
  }
  
  #--- Extract coordinates from existing ---#
  
  xy <- sf::st_coordinates(existing)
  
  vals <- terra::extract(mraster, xy)
  
  #--- extract other attributes from sampling and remove geometry attribute ---#
  
  samp_mets <- as.data.frame(existing)
  
  samp_mets <- dplyr::select(samp_mets, -geometry)
  
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
    
    sf::st_crs(samples) <- terra::crs(mraster, proj = TRUE)
    
    if (!is.null(filename)) {
      if (!is.logical(overwrite)) {
        stop("'overwrite' must be either TRUE or FALSE")
      }
      
      if (file.exists(filename) & isFALSE(overwrite)) {
        stop(glue::glue("{filename} already exists and overwrite = FALSE"))
      }
      
      sf::st_write(samples, filename, delete_layer = overwrite)
    }
    
    #--- return sf object ---#
    return(samples)
  }
}
