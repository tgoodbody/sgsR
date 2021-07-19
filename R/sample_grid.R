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
#' @param ... Additional arguments for \code{\link[sf]{st_make_grid}}.
#'
#' @return An sf object with sampled points at intersections of fishnet grid.
#'
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
                        ...) {

  #--- Set global vars ---#

  ext <- geometry <- x <-  NULL

  if (!inherits(raster, "SpatRaster")) {
    stop("'raster' must be type SpatRaster", call. = FALSE)
  }

  if (!is.numeric(gridsize)) {
    stop("'gridsize' must be type numeric")
  }

  if (gridsize < 0) {
    stop("'gridsize' must be > 0")
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical")
  }

  #--- determine crs of input raster ---#
  crs <- crs(raster, proj = TRUE)

  #--- set mraster for plotting who area in case of masking ---#

  rasterP <- raster

  if (!is.null(access)) {

    #--- error handling in the presence of 'access' ---#
    if (!inherits(access, "sf")) {
      stop("'access' must be an 'sf' object")
    }

    if (!inherits(sf::st_geometry(access), "sfc_MULTILINESTRING")) {
      stop("'access' geometry type must be 'sfc_MULTILINESTRING'")
    }

    access_buff <- mask_access(raster = raster, access = access, buff_inner = buff_inner, buff_outer = buff_outer)

    raster <- access_buff$rast
  }

  #--- convert raster extent into a polygon ---#

  sfObj <- sf::st_as_sf(terra::as.polygons(ext(raster), crs = terra::crs(raster)))

  #--- create grid and locate samples at intersections ---#

  gridSamp <- sf::st_as_sf(sf::st_make_grid(sfObj, gridsize, what = "corners", crs = terra::crs(raster), ...))

  #--- extract values from raster for each sample ---#

  gridSamp <- extract_metrics(mraster = raster, existing = gridSamp)

  #--- set geometry column and remove samples with NA values ---#

  st_geometry(gridSamp) <- "geometry"

  gridSamp <- gridSamp %>%
    dplyr::filter(!is.na(.)) %>%
    dplyr::select(-x)

  if (isTRUE(plot)) {

    #--- plot input raster and random samples ---#

    if (!is.null(access)) {
      suppressWarnings(terra::plot(rasterP[[1]]))
      suppressWarnings(terra::plot(access_buff$buff, add = T, border = c("gray30"), col = "gray10", alpha = 0.1))
      suppressWarnings(terra::plot(gridSamp, add = TRUE, col = "black"))
    } else {
      suppressWarnings(terra::plot(rasterP[[1]]))
      suppressWarnings(terra::plot(gridSamp, add = TRUE, col = "black"))
    }
  }

  if (!is.null(filename)) {
    if (!is.logical(overwrite)) {
      stop("'overwrite' must be either TRUE or FALSE")
    }

    if (file.exists(filename) & isFALSE(overwrite)) {
      stop(paste0(filename, " already exists and overwrite = FALSE"))
    }

    sf::st_write(gridSamp, filename, delete_layer = overwrite)
  }

  #--- output ---#

  return(gridSamp)
}
