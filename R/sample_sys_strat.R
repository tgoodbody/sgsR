#' Systematic stratified sampling
#'
#' @description Systematic stratified sampling with random start point and translation within a square or hexagonal tessellation for each stratum.
#'
#' @family sample functions
#' @inheritParams sample_systematic
#'
#' @param sraster spatRaster. Stratified raster with name \code{"strata"}.
#'
#' @importFrom graphics par
#'
#' @return An sf object with sampled points over unique tessellations.
#'
#' @note Specifying \code{location = "random"} can result in tessellations with no samples.
#' This results from \code{raster} have \code{NA} values at the random location chosen.
#' Using \code{force = TRUE} removes areas of \code{NA} from sampling entirely, but
#' considerably slows processing speeds. Thanks to R. Hijmans for help in debugging and
#' providing suggestions for this script.
#'
#' All stratum are sampled using random tessellation start points and translations.
#'
#' @examples
#' #--- Load raster and access files ---#
#' r <- system.file("extdata", "sraster.tif", package = "sgsR")
#' sr <- terra::rast(r)
#'
#' #--- perform grid sampling ---#
#' sample_sys_strat(
#'   sraster = sr,
#'   cellsize = 1000
#' )
#'
#' sample_sys_strat(
#'   sraster = sr,
#'   cellsize = 1000,
#'   square = FALSE,
#'   location = "corners"
#' )
#'
#' @author Tristan R.H. Goodbody
#'
#' @export

sample_sys_strat <- function(sraster,
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
                             details = FALSE) {
  strata <- NULL

  if (!inherits(sraster, "SpatRaster")) {
    stop("'sraster' must be type SpatRaster.", call. = FALSE)
  }

  if (!any("strata" %in% names(sraster))) {
    stop("'sraster' must have a layer named `strata`.", call. = FALSE)
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

  if (!is.logical(details)) {
    stop("'details' must be type logical.", call. = FALSE)
  }

  if (!is.character(location)) {
    stop("'location' must be type character.", call. = FALSE)
  } else {
    if (!any(c("centers", "corners", "random") %in% location)) {
      stop("'location' must be one of 'centers', 'corners', or 'random'.", call. = FALSE)
    }
  }

  #--- set mraster for plotting who area in case of masking ---#

  r <- sraster

  if (!missing(access)) {
    access_buff <- mask_access(raster = sraster, access = access, buff_inner = buff_inner, buff_outer = buff_outer)

    sraster <- access_buff$rast
  }

  uniquevals <- unique(terra::values(sraster, dataframe = TRUE)) %>%
    na.omit() %>%
    dplyr::arrange(strata)

  un <- nrow(uniquevals)

  tmp <- sf::st_sfc()
  class(tmp)[1] <- "sfc_POINT" # for points
  samples <- sf::st_sf(strata = integer(0), geometry = tmp, crs = terra::crs(sraster))

  tessellation <- list()


  for (i in 1:un) {
    s <- as.numeric(uniquevals[i, ])

    message(paste0("Processing strata : ", s))

    strata_m <- terra::mask(sraster,
      mask = sraster,
      maskvalues = uniquevals[i, ],
      inverse = TRUE
    )

    s <- sample_systematic(
      raster = strata_m,
      cellsize = cellsize,
      square = square,
      location = location,
      force = force,
      access = NULL,
      buff_inner = NULL,
      buff_outer = NULL,
      plot = FALSE,
      filename = NULL,
      overwrite = FALSE,
      details = TRUE
    )

    #--- get stratum values ---#
    stmp <- extract_strata(sraster = sraster, existing = s$samples, quiet = FALSE)

    samples <- rbind(samples, stmp)

    #--- add tesselation for each stratum to a list object ---#

    tessellation[[i]] <- s$tessellation %>%
      dplyr::mutate(strata = uniquevals[i, ])
  }

  if (isTRUE(plot)) {
    rasterext <- sf::st_as_sf(terra::as.polygons(terra::ext(sraster[[1]]), crs = terra::crs(sraster)))

    # number of plots
    nplots <- un

    # calculate the number of rows and columns
    ncols <- ceiling(sqrt(nplots))
    nrows <- ceiling(nplots / ncols)

    # set the mfrow parameter
    par(mfrow = c(nrows, ncols), mai = c(1, 0.1, 0.1, 0.1))

    m <- c(3.1, 3.1, 1.1, 1.1)

    # generate plots
    for (i in 1:nplots) {
      gridR <- sf::st_intersection(sf::st_geometry(tessellation[[i]]), sf::st_geometry(rasterext)) %>%
        sf::st_as_sf()

      ss <- samples %>% dplyr::filter(strata == i)

      if (!is.null(access)) {
        terra::plot(sraster, main = paste0("strata_", i), mar = m, legend = FALSE)
        terra::plot(gridR, add = TRUE, col = "transparent", border = c("blueviolet"), alpha = 0.01)
        terra::plot(access_buff$buff, add = TRUE, border = c("gray30"), col = "gray10", alpha = 0.1)
        terra::plot(ss, add = TRUE, col = "black")
      } else {
        terra::plot(r, main = paste0("strata_", i), mar = m, legend = FALSE)
        terra::plot(gridR, add = TRUE, col = "transparent", border = c("blueviolet"), alpha = 0.01)
        terra::plot(ss, add = TRUE, col = "black")
      }
    }
  }

  if (!is.null(filename)) {
    if (!is.character(filename)) {
      stop("'filename' must be a file path character string.", call. = FALSE)
    }

    if (!is.logical(overwrite)) {
      stop("'overwrite' must be type logical.", call. = FALSE)
    }

    if (file.exists(filename) & isFALSE(overwrite)) {
      stop(paste0("'", filename, "' already exists and overwrite = FALSE."), call. = FALSE)
    }

    sf::st_write(samples, filename, delete_layer = overwrite)
    message("Output samples written to disc.")
  }

  if (isTRUE(details)) {
    #--- output metrics details along with stratification raster ---#

    output <- list(samples = samples, tessellation = tessellation)

    #--- output samples dataframe ---#
    return(output)
  } else {
    #--- just output raster ---#

    return(samples)
  }
}
