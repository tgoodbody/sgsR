#' Balanced sampling
#'
#' @description Balanced raster sampling using \code{\link[BalancedSampling]{lcube}} and
#' \code{\link[SamplingBigData]{lpm2_kdtree}} methods
#'
#' @family sample functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams sample_srs
#' @param algorithm Character. One of \code{lpm2_kdtree, lcube, lcubestratified}.
#' @param p Numeric. Vector with length equal to the number of cells in \code{mraster} representing
#' the inclusion probability for each candidate sample. Default = \code{nSamp / N}, where \code{N}
#' is the number of cells.
#' @param filename Character. Path to write output samples.
#'
#' @return An sf object with \code{nSamp} samples.
#'
#' @examples
#' #--- Load raster and existing plots---#
#' r <- system.file("extdata", "mraster.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' sample_balanced(
#'   mraster = mr,
#'   nSamp = 200
#' )
#'
#' @references
#'
#' Anton Grafström and Jonathan Lisic (2019). BalancedSampling: Balanced and Spatially
#' Balanced Sampling. R package version 1.5.5. https://CRAN.R-project.org/package=BalancedSampling
#'
#' Jonathan Lisic and Anton Grafström (2018). SamplingBigData: Sampling Methods for
#' Big Data. R package version 1.0.0. https://CRAN.R-project.org/package=SamplingBigData
#'
#' Grafström, A. Lisic, J (2018). BalancedSampling: Balanced and Spatially Balanced Sampling.
#'  R package version 1.5.4. http://www.antongrafstrom.se/balancedsampling
#'
#' @author Tristan R.H. Goodbody
#'
#' @export

sample_balanced <- function(mraster,
                            nSamp,
                            algorithm = "lpm2_kdtree",
                            p = NULL,
                            access = NULL,
                            buff_inner = NULL,
                            buff_outer = NULL,
                            plot = FALSE,
                            filename = NULL,
                            overwrite = FALSE) {
  #--- Set global vars ---#
  x <- y <- X <- Y <- strata <- NULL

  #--- Error management ---#
  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster.", call. = FALSE)
  }

  if (!is.numeric(nSamp)) {
    stop("'nSamp' must be type numeric.", call. = FALSE)
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical.", call. = FALSE)
  }

  if (!is.character(algorithm)) {
    stop("'algorithm' must be type character.", call. = FALSE)
  }

  #--- list all available algorithms to determine if a valid one has been supplied ---#
  algs <- c("lpm2_kdtree", "lcube", "lcubestratified")

  if (!algorithm %in% algs) {
    stop("Unknown algorithm specified. Please use one of 'lpm2_kdtree', 'lcube', 'lcubestratified'.", call. = FALSE)
  }

  #--- determine crs of input mraster ---#
  crs <- terra::crs(mraster, proj = TRUE)

  #--- set mraster for plotting who area in case of masking ---#

  mrasterP <- mraster

  if (!is.null(access)) {
    #--- buffer roads and mask ---#

    access_buff <- mask_access(
      raster = mraster,
      access = access,
      buff_inner = buff_inner,
      buff_outer = buff_outer
    )

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
      stop("'p' must be type numeric.", call. = FALSE)
    }

    if (length(p) != N) {
      stop(paste0("'p' have a length of ", N, "."), call. = FALSE)
    }
  }


  if (algorithm == "lpm2_kdtree") {
    sampled <- SamplingBigData::lpm2_kdtree(prob = p, x = vals_m)
  }

  if (algorithm == "lcube") {
    sampled <- BalancedSampling::lcube(prob = p, Xspread = vals_m, Xbal = cbind(p))
  }

  if (algorithm == "lcubestratified") {
    if (!"strata" %in% names(mraster)) {
      stop("'mraster' must have a variable named 'strata' to use the 'lcubestratified' algorithm.", call. = FALSE)
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

  }

  samples <- vals[sampled, ]

  #--- convert coordinates to a spatial points object ---#
  samples <- dplyr::select(samples, X, Y) %>%
    as.data.frame() %>%
    sf::st_as_sf(., coords = c("X", "Y"), crs = crs)

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

  #--- write outputs if desired ---#
  write_samples(samples = samples, filename = filename, overwrite = overwrite)

  #--- output samples sf ---#

  return(samples)
}
