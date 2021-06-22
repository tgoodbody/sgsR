#' Conditioned Latin Hypercube Sampling
#' @family sample functions
#'
#' @importFrom magrittr %>%
#' @importFrom methods is
#' 
#' @return An sf object with \code{n} stratified samples.
#' 
#' @export

sample_clhs <- function(mraster = mraster,
                        n = NULL,
                        access = NULL,
                        cost.access = 
                        plot = FALSE) 
{
  
  #--- Error management ---#
  
  #--- Error management ---#
  if (!inherits(mraster, "SpatRaster"))
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  
  if (!is.numeric(n))
    stop("'n' must be type numeric")
  
  if (!is.logical(plot))
    stop("'plot' must be type logical")
  
  
  if (is.na(crs(mraster)))
    stop("'mraster' does not have a coordinate system")
  
  #--- determine crs of input mraster ---#
  crs <- crs(mraster)
  
  if (!is.null(access)) {
  
    #--- error handling in the presence of 'access' ---#
    if (!inherits(access,"sf"))
      stop("'access' must be an 'sf' object")
    
    if(!inherits(sf::st_geometry(access),"sfc_MULTILINESTRING"))
      stop("'access' geometry type must be 'sfc_MULTILINESTRING'")
    
    #--- list all buffers to catch NULL values within error handling ---#
    buffers <- list(buff_inner, buff_outer)
    
    #--- error handling in the presence of 'access' ---#
    if (any(vapply(buffers, is.null, TRUE)))
      stop("All 'buff_*' paramaters must be provided when 'access' is defined.")
    
    if (!any(vapply(buffers, is.numeric, FALSE)))
      stop("All 'buff_*' paramaters must be type numeric")
    
      message(
        paste0(
          "An access layer has been provided. An internal buffer of ",
          buff_inner,
          " m and an external buffer of ",
          buff_outer,
          " m have been applied"
        )
      )
      


  }
  
  #--- convert vectors to spatVector to synergize with terra raster functions---#
  roads <- terra::vect(roads)
  
  #--- make access buffer with user defined values ---#
  
  buff_in <- terra::buffer(x = roads,
                           width = 50)
  
  #--- make difference and aggregate inner and outer buffers to prevent sampling too close to access ---#
  
  t <- terra::rasterize(buff_in,mraster[[1]])
  
  mmm <- mask(x = mraster[[1]], mask = t, inverse = TRUE)
  
  
  
  #--- extract XY coordinates from raster ---#
  vals <- terra::as.data.frame(wall_poly[[1:3]],xy=TRUE) %>%
    dplyr::rename(X = x,
                  Y = y)
  
  
  roads <- terra::vect(roads)
  
  buff_in <- terra::buffer(x = roads,
                           width = 100)
  
  
  wall_poly <- terra::mask(wall_poly[[1]], mask = buff_in, inverse = TRUE)

  
  kk <- terra::distance(mraster[[1]],roads)
  
  kk.m <- mask(kk,wall_poly[[1]])
  
  terra::as.data.frame(kk.m,xy=TRUE) %>%
    dplyr::rename(X = x,
                  Y = y)
  
  

  
}
