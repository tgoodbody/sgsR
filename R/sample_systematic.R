#' Systematic sampling
#'
#' @description Systematic sampling within a square or hexagonal tessellation.
#'
#' @family sample functions
#'
#' @param raster spatRaster. Raster used to define extent of fishnet grid.
#' @param cellsize Numeric. Desired cellsize for tessellation.
#' @param square Logical. Tessellation shape. Default is regular square grid,
#' if \code{FALSE} hexagons are returned.
#' @param centers Logical. Sample location within tessellation. Default (\code{TRUE})
#' returns samples at centers. If \code{FALSE}, corners are returned.
#' @param access sf. Road access network - must be lines.
#' @param buff_inner Numeric. Inner buffer boundary specifying distance
#'  from access where plots cannot be sampled.
#' @param buff_outer Numeric. Outer buffer boundary specifying distance
#'  from access where plots can be sampled.
#' @param plot Logical. Plots output strata raster with samples.
#' @param filename Character. Path to write output samples.
#' @param overwrite Logical. Choice to overwrite existing \code{filename} if it exists.
#' @param details Logical. If \code{FALSE} (default) output is sf object of
#' systematic samples. If \code{TRUE} returns a list of sf objects where \code{tessellation}
#' is the tessellation grid for sampling, and \code{samples} are the systematic samples.
#' @param ... Additional arguments for \code{\link[sf]{st_make_grid}}.
#'
#' @return An sf object with sampled points over a tessellation.
#'
#' @examples
#' #--- Load raster and access files ---#
#' r <- system.file("extdata", "kmeans.tif", package = "sgsR")
#' sr <- terra::rast(r)
#'
#' a <- system.file("extdata", "roads.shp", package = "sgsR")
#' ac <- sf::st_read(a)
#'
#' e <- system.file("extdata", "existing.shp", package = "sgsR")
#' e <- sf::st_read(e)
#'
#' #--- perform grid sampling ---#
#' sample_systematic(raster = sr, 
#'             cellsize = 200, 
#'             plot = TRUE)
#'             
#' sample_systematic(raster = sr, 
#'             cellsize = 200,
#'             square = FALSE,
#'             centers = FALSE,
#'             plot = TRUE)
#' 
#' sample_systematic(raster = sr, 
#'             cellsize = 100, 
#'             access = ac,
#'             buff_inner = 50,
#'             buff_outer = 200)
#' 
#' sample_systematic(raster = sr,
#'             cellsize = 200,
#'             access = ac,
#'             buff_inner = 100,
#'             buff_outer = 400,
#'             filename = tempfile(fileext = ".shp"),
#'             plot = TRUE)
#'             
#' @author Tristan R.H. Goodbody
#'
#' @export


sample_systematic <- function(raster,
                              cellsize,
                              square = TRUE,
                              centers = TRUE,
                              access = NULL,
                              buff_inner = NULL,
                              buff_outer = NULL,
                              plot = FALSE,
                              filename = NULL,
                              overwrite = FALSE,
                              details = FALSE,
                              ...) {
  
  #--- Set global vars ---#
  
  ext <- geometry <- x <- NULL
  
  if (!inherits(raster, "SpatRaster")) {
    stop("'raster' must be type SpatRaster", call. = FALSE)
  }
  
  if (!is.numeric(cellsize)) {
    stop("'cellsize' must be type numeric")
  }
  
  if (cellsize < 0) {
    stop("'cellsize' must be > 0")
  }
  
  if (!is.logical(plot)) {
    stop("'plot' must be type logical")
  }
  
  if (!is.logical(square)) {
    stop("'square' must be type logical")
  }
  
  if (!is.logical(centers)) {
    stop("'centers' must be type logical")
    
  } else {
    
    #--- depending on whether 'centers' is true or false allocate choice to 'location' ---#
    
    if(isTRUE(centers)){
      location <-  "centers"
    } else {
      location <- "corners"
    }
    
  }
  
  #--- determine crs of input raster ---#
  crs <- terra::crs(raster, proj = TRUE)
  
  #--- set mraster for plotting who area in case of masking ---#
  
  rasterP <- raster
  
  if (!is.null(access)) {
    
    #--- error handling in the presence of 'access' ---#
    if (!inherits(access, "sf")) {
      stop("'access' must be an 'sf' object")
    }
    
    if (!inherits(sf::st_geometry(access), "sfc_MULTILINESTRING") && !inherits(sf::st_geometry(access), "sfc_LINESTRING")) {
      stop("'access' geometry type must be 'LINESTRING' or 'MULTILINESTRING'")
    }
    
    access_buff <- mask_access(raster = raster, access = access, buff_inner = buff_inner, buff_outer = buff_outer)
    
    raster <- access_buff$rast
  }
  
  #--- convert raster extent into a polygon ---#
  
  sfObj <- sf::st_as_sf(terra::as.polygons(terra::ext(raster), crs = terra::crs(raster)))
  
  #--- create grid and locate samples ---#
  
  samples <- sf::st_as_sf(sf::st_make_grid(sfObj, cellsize, square = square, what = location, crs = terra::crs(raster), ...)) %>%
    dplyr::rename(geometry = x) %>%
    #--- need to extract a metric to determine if values are NA ---#
    extract_metrics(mraster = raster[[1]], existing = .) %>%
    #--- remove samples with NA ---#
    dplyr::filter(!is.na(.)) %>%
    dplyr::select(geometry)
  
  #--- create tessellation ---#
  
  grid <- sf::st_as_sf(sf::st_make_grid(sfObj, cellsize, square = square, what = "polygons", crs = terra::crs(raster), ...))
  

  if (isTRUE(plot)) {
    
    #--- plot input raster and random samples ---#
    
    if (!is.null(access)) {
      suppressWarnings(terra::plot(rasterP[[1]]))
      suppressWarnings(terra::plot(grid, add = TRUE, border = c("blueviolet"), alpha = 0.01))
      suppressWarnings(terra::plot(access_buff$buff, add = T, border = c("gray30"), col = "gray10", alpha = 0.1))
      suppressWarnings(terra::plot(samples, add = TRUE, col = "black"))
    } else {
      suppressWarnings(terra::plot(rasterP[[1]]))
      suppressWarnings(terra::plot(grid, add = TRUE, border = c("blueviolet"), alpha = 0.01))
      suppressWarnings(terra::plot(samples, add = TRUE, col = "black"))
    }
  }
  
  if (!is.null(filename)) {
    if (!is.logical(overwrite)) {
      stop("'overwrite' must be either TRUE or FALSE")
    }
    
    if (file.exists(filename) & isFALSE(overwrite)) {
      stop(glue::glue('{filename} already exists and overwrite = FALSE'))
    }
    
    sf::st_write(samples, filename, delete_layer = overwrite)
  }
  
  if (isTRUE(details)) {
    
    #--- output metrics details along with stratification raster ---#
    
    output <- list(samples = samples, tessellation = grid)
    
    #--- output samples dataframe ---#
    return(output)
  } else {
    
    #--- just output raster ---#
    
    return(samples)
  }
}
