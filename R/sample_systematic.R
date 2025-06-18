#' Systematic sampling
#'
#' @description Systematic sampling with random start point and translation within a square or hexagonal tessellation.
#'
#' @family sample functions
#'
#' @param raster spatRaster. Raster used to define extent of fishnet grid.
#' @param cellsize Numeric. Desired cellsize for tessellation.
#' @param square Logical. Tessellation shape. Default is regular square grid,
#' if \code{FALSE} hexagons are used.
#' @param location Character. Sample location within tessellation. \code{Default = "centers"})
#' returns samples at tessellation centers, \code{"corners"} - corners of tessellation are returned,
#' \code{"random"} - samples are randomly located within tessellations.
#' @param force Logical. Only applies when \code{location = "random"}. If \code{TRUE}, random samples are
#' forced to fall in areas where \code{raster} does not have \code{NA} values. This will considerably slow processing.
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
#' @param ... Additional arguments for \code{\link[sf]{st_make_grid}}. Options include \code{offset}
#' to offset grid by providing lower left coordinates.
#'
#' @return An sf object with sampled points over a tessellation.
#'
#' @note Specifying \code{location = "random"} can result in tessellations with no samples.
#' This results from \code{raster} have \code{NA} values at the random location chosen.
#' Using \code{force = TRUE} removes areas of \code{NA} from sampling entirely, but
#' considerably slows processing speeds. Thanks to R. Hijmans for help in debugging and
#' providing suggestions for this script.
#'
#' @examples
#' #--- Load raster and access files ---#
#' r <- system.file("extdata", "mraster.tif", package = "sgsR")
#' mr <- terra::rast(r)
#' #--- perform grid sampling ---#
#' sample_systematic(
#'   raster = mr,
#'   cellsize = 1000
#' )
#'
#' sample_systematic(
#'   raster = mr,
#'   cellsize = 1000,
#'   location = "corners",
#'   plot = TRUE
#' )
#'
#' sample_systematic(
#'   raster = mr,
#'   cellsize = 1000,
#'   square = FALSE,
#'   location = "random",
#'   plot = TRUE
#' )
#'
#' @author Tristan R.H. Goodbody, Lukas Winiwarter
#'
#' @export

sample_systematic <- function(raster,
                              cellsize,
                              square = TRUE,
                              location = "centers",
                              force = FALSE,
                              access = NULL,
                              buff_inner = NULL,
                              buff_outer = NULL,
                              plot = FALSE,
                              filename = NULL,
                              overwrite = FALSE,
                              details = FALSE,
                              ...) {
  #--- Set global vars ---#

  ext <- geometry <- x <- pass <- overlap <- NULL

  if (!inherits(raster, "SpatRaster")) {
    stop("'raster' must be type SpatRaster.", call. = FALSE)
  }

  if (!is.numeric(cellsize)) {
    stop("'cellsize' must be type numeric.", call. = FALSE)
  }

  if (cellsize < 0) {
    stop("'cellsize' must be > 0.", call. = FALSE)
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical.", call. = FALSE)
  }

  if (!is.logical(square)) {
    stop("'square' must be type logical.", call. = FALSE)
  }

  if (!is.character(location)) {
    stop("'location' must be type character.", call. = FALSE)
  } else {
    if (!any(c("centers", "corners", "random") %in% location)) {
      stop("'location' must be one of 'centers', 'corners', or 'random'.", call. = FALSE)
    }
  }

  #--- determine crs of input raster ---#
  crs <- terra::crs(raster)

  #--- set mraster for plotting who area in case of masking ---#

  r <- raster

  if (!is.null(access)) {
    access_buff <- mask_access(raster = raster, access = access, buff_inner = buff_inner, buff_outer = buff_outer)

    raster <- access_buff$rast
  }

  #--- convert raster extent into a polygon ---#

  res <- terra::xres(raster)
  e <- as.vector(terra::ext(raster))

  #--- add randomness to grid lower left coordinate locations ---#

  xminc <- e[1] + (res * sample(-100:0, 1))
  yminc <- e[3] + (res * sample(-100:0, 1))

  #--- generate raster extent polygon ---#

  rasterext <- sf::st_as_sf(terra::as.polygons(terra::ext(raster[[1]]), crs = crs))

  #--- add random LL corner to raster extent to mark start point of systematic grid ---#
  pol <- terra::as.polygons(terra::ext(xminc, e[2], yminc, e[4]), crs = crs)

  #--- generate a seperate sf component ---#
  polsf <- sf::st_as_sf(x = pol, crs = crs)

  #--- random translation value ---#

  randRot <- runif(1, 0, 360)

  #--- spin pol and get extent for grid making to cover entire raster ---#

  se <- terra::spin(pol, randRot) %>% terra::ext()

  #--- create grid based on outer extents of randomly translated raster extent ---#

  grid <- terra::vect(sf::st_make_grid(
    x = se,
    cellsize = cellsize,
    square = square,
    what = "polygons",
    crs = crs,
    ...
  ))

  #--- randomly spin the grid ---#
  se <- terra::spin(grid, randRot)

  terra::crs(se, warn = FALSE) <- crs

  #--- check which polygons are fully overlapping with the raster extent ---#

  gridR <- sf::st_as_sf(se) %>%
    dplyr::mutate(overlap = lengths(sf::st_intersects(., rasterext))) %>%
    dplyr::filter(overlap == 1)


  #--- check to make sure that samples intersect raster extent (cellsize check) ---#

  if (is.null(nrow(gridR)) || nrow(gridR) == 0) {
    stop("No samples intersect with 'raster'. Ensure 'cellsize' makes sense.", call. = FALSE)
  }

  if (isTRUE(force)) {
    if (location != "random") {
      stop("'location' must be 'random' when 'force = TRUE'", call. = FALSE)
    }

    message("Forcing samples to fall in non NA locations.")

    #--- force "random" samples to not fall in areas of no data ---#

    #--- generate NA raster mask ---#
    m <- terra::setValues(raster[[1]], NA)
    m[is.na(raster[[1]])] <- 1

    #--- convert non-NA areas to polygon ---#
    p <- sf::st_as_sf(terra::as.polygons(m)) %>%
      dplyr::select(geometry) %>%
      sf::st_combine()

    #--- if there is overlapping or slivers fix the issues ---#
    if (isFALSE(sf::st_is_valid(p))) {
      p <- p %>%
        sf::st_make_valid()
    }

    p_diff <- sf::st_difference(sf::st_geometry(gridR), sf::st_geometry(p)) %>%
      sf::st_as_sf() %>%
      sf::st_intersection(., rasterext) %>%
      dplyr::filter(sf::st_is(., c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION")))

    #--- sampling --#
    samples <- sf::st_sample(p_diff, size = c(1, 1), type = "random") %>%
      sf::st_sf() %>%
      dplyr::mutate(overlap = lengths(sf::st_intersects(., rasterext))) %>%
      dplyr::filter(overlap == 1)

    #--- check to make sure that samples intersect raster extent (cellsize check) ---#

    if (is.null(nrow(samples)) || nrow(samples) == 0) {
      stop("No samples intersect with 'raster'. Ensure 'cellsize' makes sense.", call. = FALSE)
    }

    samples <- samples %>%
      extract_metrics(mraster = raster[[1]], existing = ., quiet = TRUE) %>%
      stats::na.omit() %>%
      dplyr::select(geometry)
  } else {
    #--- random sampling within tesselations ---#

    if (location == "random") {
      location <- "centers"

      #--- maximum number of samples that can be selected ---#

      gridn <- nrow(gridR)

      #--- determine maximum distance sample can be moved in X / Y to remain within each tessellation ---#

      radius <- cellsize / 2

      if (square == TRUE) {
        vals <- data.frame(
          X = runif(gridn, -radius, radius),
          Y = runif(gridn, -radius, radius)
        )
      } else {
        #--- parameters for hexagon random sampling ---#
        c30 <- sqrt(3) / 2

        tests <- gridn * 100

        X <- runif(tests, -1, 1)
        Y <- runif(tests, -c30, c30)

        xy <- data.frame(X = X, Y = Y)

        #--- test whether the sample will be within the bounds of the hexagon ---#
        xy$pass <- abs(xy$X) < 1 - .5 * (abs(xy$Y) / c30)

        #--- filter only values with TRUE in $pass ---#
        vals <- xy %>%
          dplyr::filter(pass == TRUE) %>%
          dplyr::slice_sample(., n = nrow(gridR)) %>%
          dplyr::mutate(
            X = X * radius,
            Y = Y * radius
          )
      }

      #--- create grid and locate samples ---#

      samples <- sf::st_centroid(sf::st_geometry(gridR)) %>%
        sf::st_as_sf() %>%
        dplyr::rename(geometry = x) %>%
        sf::st_coordinates(.) %>%
        as.data.frame() %>%
        #--- apply random movement by row ---#
        dplyr::mutate(
          X = (vals$X * cos(randRot)) + (vals$Y * sin(randRot)) + X,
          Y = (vals$X * -sin(randRot)) + (vals$Y * cos(randRot)) + Y
        ) %>%
        sf::st_as_sf(.,
          coords = c("X", "Y"),
          crs = crs
        ) %>%
        dplyr::mutate(overlap = lengths(sf::st_intersects(., rasterext))) %>%
        dplyr::filter(overlap == 1)

      #--- check to make sure that samples intersect raster extent (cellsize check) ---#

      if (is.null(nrow(samples)) || nrow(samples) == 0) {
        stop("No samples intersect with 'raster'. Ensure 'cellsize' makes sense.", call. = FALSE)
      }

      samples <- samples %>%
        extract_metrics(mraster = raster[[1]], existing = ., quiet = TRUE) %>%
        stats::na.omit() %>%
        dplyr::select(geometry)
    } else if (location == "corners") {
      samples <- sf::st_cast(sf::st_geometry(gridR), "POINT") %>%
        sf::st_as_sf() %>%
        dplyr::rename(geometry = x) %>%
        dplyr::mutate(overlap = lengths(sf::st_intersects(., rasterext))) %>%
        dplyr::filter(overlap == 1)

      #--- check to make sure that samples intersect raster extent (cellsize check) ---#

      if (is.null(nrow(samples)) || nrow(samples) == 0) {
        stop("No samples intersect with 'raster'. Ensure 'cellsize' makes sense.", call. = FALSE)
      }

      samples <- samples %>%
        extract_metrics(mraster = raster[[1]], existing = ., quiet = TRUE) %>%
        stats::na.omit() %>%
        dplyr::select(geometry)
    } else if (location == "centers") {
      samples <- sf::st_centroid(sf::st_geometry(gridR)) %>%
        sf::st_as_sf() %>%
        dplyr::rename(geometry = x) %>%
        dplyr::mutate(overlap = lengths(sf::st_intersects(., rasterext))) %>%
        dplyr::filter(overlap == 1)

      #--- check to make sure that samples intersect raster extent (cellsize check) ---#

      if (is.null(nrow(samples)) || nrow(samples) == 0) {
        stop("No samples intersect with 'raster'. Ensure 'cellsize' makes sense.", call. = FALSE)
      }

      samples <- samples %>%
        extract_metrics(mraster = raster[[1]], existing = ., quiet = TRUE) %>%
        stats::na.omit() %>%
        dplyr::select(geometry)
    }
  }

  if (isTRUE(plot)) {
    #--- plot input raster and random samples ---#

    gridR <- sf::st_intersection(sf::st_geometry(gridR), sf::st_geometry(rasterext)) %>%
      sf::st_as_sf()

    if (!is.null(access)) {
      terra::plot(r[[1]])
      terra::plot(gridR, add = TRUE, col = "transparent", border = c("blueviolet"), alpha = 0.01)
      terra::plot(access_buff$buff, add = T, border = c("gray30"), col = "gray10", alpha = 0.1)
      terra::plot(samples, add = TRUE, col = "black")
    } else {
      terra::plot(r[[1]])
      terra::plot(gridR, add = TRUE, col = "transparent", border = c("blueviolet"), alpha = 0.01)
      terra::plot(samples, add = TRUE, col = "black")
    }
  }

  #--- write outputs if desired ---#
  write_samples(samples = samples, filename = filename, overwrite = overwrite)

  if (isTRUE(details)) {
    #--- output metrics details along with stratification raster ---#

    output <- list(samples = samples, tessellation = gridR)

    #--- output samples dataframe ---#
    return(output)
  } else {
    #--- just output raster ---#

    return(samples)
  }
}
