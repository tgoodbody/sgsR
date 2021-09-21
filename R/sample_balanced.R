#' Balanced sampling
#'
#' @description Balanced raster sampling using \code{\link[BalancedSampling]{lpm2}} and
#' \code{\link[SamplingBigData]{lpm2_kdtree}} methods
#'
#' @family sample functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams sample_srs
#' @param algorithm Character. One of \code{lpm2 lcube lcubestratified}
#' @param p Numeric. Vector with length equal to the number of cells in \code{mraster} representing
#' the inclusion probability for each candidate sample. Default = \code{nSamp / N}, where \code{N}
#' is the number of cells.
#'
#' @return An sf object with \code{nSamp} randomly sampled points.
#'
#' @examples
#' \dontrun{
#' #--- Load raster and existing plots---#
#' r <- system.file("extdata", "wall_metrics_small.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' a <- system.file("extdata", "roads.shp", package = "sgsR")
#' ac <- sf::st_read(a)
#'
#' sample_balanced(
#'   mraster = mr,
#'   nSamp = 200,
#'   plot = TRUE
#' )
#'
#' sample_balanced(
#'   mraster = mr,
#'   nSamp = 100,
#'   algorithm = "lcube",
#'   access = ac,
#'   buff_inner = 50,
#'   buff_outer = 200
#' )
#' }
#'
#' @references
#'
#' Anton Grafstrom and Jonathan Lisic (2019). BalancedSampling: Balanced and Spatially
#' Balanced Sampling. R package version 1.5.5. https://CRAN.R-project.org/package=BalancedSampling
#'
#' Jonathan Lisic and Anton Grafstrom (2018). SamplingBigData: Sampling Methods for
#' Big Data. R package version 1.0.0. https://CRAN.R-project.org/package=SamplingBigData
#'
#' @author Tristan R.H. Goodbody
#'
#' @export

sample_balanced <- function(mraster,
                            nSamp,
                            algorithm = "lpm2",
                            p = NULL,
                            access = NULL,
                            buff_inner = NULL,
                            buff_outer = NULL,
                            plot = FALSE,
                            filename = NULL,
                            overwrite = FALSE) {

  #--- check for required packages ---#
  if (!requireNamespace("BalancedSampling", quietly = TRUE)) {
    stop("Package \"BalancedSampling\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  #--- Set global vars ---#
  x <- y <- X <- Y <- strata <- NULL

  #--- Error management ---#
  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  }

  if (!is.numeric(nSamp)) {
    stop("'nSamp' must be type numeric")
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical")
  }

  if (!is.character(algorithm)) {
    stop("'algorith' must be type character")
  }

  #--- list all available algorithms to determine if a valid one has been supplied ---#
  algs <- c("lpm2", "lcube", "lcubestratified")

  if (!algorithm %in% algs) {
    stop("Unknown algorithm specified. Please use one of 'lpm2' 'lcube' 'lcubestratified'")
  }

  ######################################
  ## DETERMINE NULL / NA SYNTAX FOR CRS##
  ######################################

  if (is.na(terra::crs(mraster, proj = TRUE))) {
    stop("'mraster' does not have a coordinate system")
  }

  #--- determine crs of input mraster ---#
  crs <- terra::crs(mraster, proj = TRUE)

  #--- set mraster for plotting who area in case of masking ---#

  mrasterP <- mraster

  if (!is.null(access)) {

    #--- error handling in the presence of 'access' ---#
    if (!inherits(access, "sf")) {
      stop("'access' must be an 'sf' object")
    }

    if (!inherits(sf::st_geometry(access), "sfc_MULTILINESTRING") && !inherits(sf::st_geometry(access), "sfc_LINESTRING")) {
      stop("'access' geometry type must be 'LINESTRING' or 'MULTILINESTRING'")
    }

    #--- buffer roads and mask ---#

    access_buff <- mask_access(raster = mraster, access = access, buff_inner = buff_inner, buff_outer = buff_outer)

    mraster <- access_buff$rast
  }

  #--- extract XY coordinates from raster ---#
  vals <- terra::as.data.frame(mraster, xy = TRUE) %>%
    dplyr::rename(
      X = x,
      Y = y
    )

  #--- # drop x,y matrix of auxiliary variables ---#
  vals_m <- as.matrix(dplyr::select(vals, -X, -Y))

  N <- nrow(vals)

  #--- inclusion probability ---#

  if (is.null(p)) {

    #--- if 'p' is not defined use the default ---#

    p <- rep(nSamp / N, N)
  } else {
    if (!is.numeric(p)) {
      stop("'p' must be type numeric")
    }
  }


  if (algorithm == "lpm2") {

    #--- check for required packages ---#

    if (!requireNamespace("SamplingBigData", quietly = TRUE)) {
      stop("Package \"SamplingBigData\" needed for the 'lpm2' algorithm. Please install it.",
        call. = FALSE
      )
    }

    sampled <- SamplingBigData::lpm2_kdtree(prob = p, x = vals_m)
  }

  if (algorithm == "lcube") {
    sampled <- BalancedSampling::lcube(prob = p, Xspread = vals_m, Xbal = cbind(p))
  }

  if (algorithm == "lcubestratified") {
    if (!"strata" %in% names(mraster)) {
      stop("'mraster' must have a variable named 'strata' to use the 'lcubestratified' algorithm")
    }

    #--- create indices for all, NA, and valid sampling candidates ---#

    strata_v <- as.vector(vals$strata)

    #--- remove strata as a sampling variable and convert to matrix ---#

    #--- # drop x,y matrix of auxiliary variables ---#
    vals_m <- as.matrix(dplyr::select(vals, -X, -Y, -strata))

    #--- generates a binary output where 0 is not sampled and 1 is sampled ---#

    sampled <- BalancedSampling::lcubestratified(
      prob = p,
      Xspread = vals_m,
      Xbal = cbind(p),
      integerStrata = strata_v
    )

    #--- extract all 1 (sampled) cells ---#

    sampled <- (1:N)[sampled == 1]
  }

  samples <- vals[sampled, ]

  #--- convert coordinates to a spatial points object ---#
  samples <- dplyr::select(samples, X, Y) %>%
    as.data.frame() %>%
    sf::st_as_sf(., coords = c("X", "Y"))

  #--- assign mraster crs to spatial points object ---#
  sf::st_crs(samples) <- crs

  if (isTRUE(plot)) {

    #--- plot input mraster and random samples ---#
    if (!is.null(access)) {
      terra::plot(mrasterP[[1]])
      suppressWarnings(terra::plot(access_buff$buff, add = T, border = c("gray30"), col = "gray10", alpha = 0.1))
      suppressWarnings(terra::plot(samples, add = T, col = "black"))
    } else {
      terra::plot(mrasterP[[1]])
      suppressWarnings(terra::plot(samples, add = T, col = "black"))
    }
  }

  if (!is.null(filename)) {
    if (!is.logical(overwrite)) {
      stop("'overwrite' must be either TRUE or FALSE")
    }

    if (file.exists(filename) & isFALSE(overwrite)) {
      stop(glue::glue("{filename} already exists and overwrite = FALSE"))
    }

    sf::st_write(samples, filename, delete_layer = overwrite)
  }

  #--- output samples sf ---#

  return(samples)
}
