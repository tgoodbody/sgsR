#' Grid sampling
#'  
#' @description Sample landscape using a fishnet grid pattern.
#' 
#' @family sample functions
#'
#' @param raster spatRaster. Raster used to define extent of fishnet grid
#' @param gridsize Numeric. Desired distance between samples
#' @param access sf. Road access network - must be lines.
#' @param buff_inner Numeric. Inner buffer boundary specifying distance
#'  from access where plots cannot be sampled.
#' @param buff_outer Numeric. Outer buffer boundary specifying distance
#'  from access where plots can be sampled.
#' @param plot Logical. Plots output strata raster with samples.
#' @param filename Character. 
#' @param filename Character. Path to write output samples.
#' @param overwrite Logical. Choice to overwrite existing \code{filename} if it exists.
#' @param ... Additional arguments for \link[sf]{st_make_grid}.
#' 
#' @return An sf object with sampled points at intersections of fishnet grid.
#' 
#' @importFrom magrittr %>%
#' @importFrom methods is
#' 
#' @export


sample_grid <- function(raster,
                        gridsize,
                        access = NULL,
                        buff_inner = NULL,
                        buff_outer = NULL,
                        plot = FALSE,
                        filename = NULL,
                        overwrite = FALSE,
                        ...)
  {
  
  if (!inherits(raster, "SpatRaster"))
    stop("'raster' must be type SpatRaster", call. = FALSE)
  
  if (!is.numeric(gridsize))
    stop("'gridsize' must be type numeric")
  
  if (gridsize < 0)
    stop("'gridsize' must be > 0")
  
  if (!is.logical(plot))
    stop("'plot' must be type logical")
  
  #--- determine crs of input raster ---#
  crs <- crs(raster)
  
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
    
    #--- convert vectors to spatVector to synergize with terra raster functions---#
    roads <- terra::vect(access)
    
    #--- make access buffer with user defined values ---#
    
    buff_in <- terra::buffer(x = roads,
                             width = buff_inner)
    
    buff_out <- terra::buffer(x = roads,
                              width = buff_outer)
    
    #--- make difference and aggregate inner and outer buffers to prevent sampling too close to access ---#
    buffer <- terra::aggregate(buff_out - buff_in)
    
    raster <- terra::mask(raster, mask = buffer)
    
  }
  
  #--- convert raster extent into a polygon ---#
  
  sfObj <- sf::st_as_sf(terra::as.polygons(ext(raster), crs = terra::crs(raster)))
  
  #--- create grid and locate samples at intersections ---#
  
  gridSamp <- sf::st_as_sf(sf::st_make_grid(sfObj, gridsize, what = "corners", crs = terra::crs(raster), ...))
  
  #--- extract values from raster for each sample ---#
  
  gridSamp <- extract_metrics(mraster = raster, samples = gridSamp)
  
  #--- remove samples with NA values ---#
  
  gridSamp <- gridSamp %>% 
    dplyr::filter(!is.na(.)) %>%
    dplyr::select(-geometry)
  
  if(isTRUE(plot)){
    
    #--- plot input raster and random samples ---#
    
    terra::plot(raster)
    terra::plot(gridSamp, add = TRUE, col = "black")
    
  }
  
  if (!is.null(filename)){
    
    if(!is.logical(overwrite))
      stop("'overwrite' must be either TRUE or FALSE")
    
    if(file.exists(filename) & isFALSE(overwrite))
      stop(paste0(filename, " already exists and overwrite = FALSE"))
    
    sf::st_write(gridSamp, filename, delete_layer = overwrite)
    
  }
  
  #--- output ---#
  
  return(gridSamp)
  
}
