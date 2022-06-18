#' k-means stratification
#'
#' @description Stratify metrics raster using \code{\link[stats]{kmeans}} algorithm
#' @family stratify functions
#'
#' @inheritParams strat_breaks
#'
#' @param mraster spatRaster. ALS metrics raster.
#' @param nStrata Character. Number of desired strata.
#' @param iter Numeric. The maximum number of iterations allowed.
#' @param algorithm Character. \code{Lloyd} (default) or
#' \code{MacQueen} algorithms.
#' @param center Logical. Value indicating whether the variables should be shifted to be zero centered.
#' @param scale Logical. Value indicating whether the variables should be scaled to have unit variance.
#' @param plot Logical. Plots output strata raster and visualized
#'  strata with boundary dividers.
#' @param details Logical. If \code{FALSE} (default) output is only
#' stratification raster. If \code{TRUE} return a list where \code{$details} is additional 
#' stratification information and \code{$raster} is the output stratification spatRaster.
#' @param ... Additional arguments to be passed to \code{\link[stats]{kmeans}} function.
#'
#' @return output stratification \code{spatRaster}, or a list when \code{details = TRUE}.
#'
#' @examples
#' #--- Load raster and access files ---#
#' r <- system.file("extdata", "mraster.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' #--- perform stratification using k-means ---#
#' kmeans <- strat_kmeans(
#'   mraster = mr,
#'   nStrata = 5
#' )
#'
#' kmeans <- strat_kmeans(
#'   mraster = mr,
#'   nStrata = 5,
#'   iter = 1000,
#'   algorithm = "MacQueen",
#'   details = TRUE
#' )
#'
#' kmeans <- strat_kmeans(
#'   mraster = mr,
#'   nStrata = 5,
#'   iter = 1000,
#'   plot = TRUE,
#'   filename = tempfile(fileext = ".tif"),
#'   overwrite = TRUE
#' )
#' @author Tristan R.H. Goodbody
#'
#' @export

strat_kmeans <- function(mraster,
                         nStrata,
                         iter = 500,
                         algorithm = "Lloyd",
                         center = TRUE,
                         scale = TRUE,
                         plot = FALSE,
                         details = FALSE,
                         filename = NULL,
                         overwrite = FALSE,
                         ...) {

  #--- Error management ---#

  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster.", call. = FALSE)
  }

  if (!is.numeric(nStrata)) {
    stop("'nStrata' must be type numeric.", call. = FALSE)
  }

  if (!is.numeric(iter)) {
    stop("'iter' must be type numeric.", call. = FALSE)
  }

  if (!is.character(algorithm)) {
    stop("'algorithm' must be type character.", call. = FALSE)
  }

  if (algorithm != "Lloyd" && algorithm != "MacQueen") {
    stop("Unknown algorithm '", algorithm, "' selected. Please use 'Lloyd' (default) or 'MacQueen'.", call. = FALSE)
  }

  if (!is.logical(center)) {
    stop("'center' must be type logical.", call. = FALSE)
  }

  if (!is.logical(scale)) {
    stop("'scale' must be type logical.", call. = FALSE)
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical.", call. = FALSE)
  }

  if (!is.logical(details)) {
    stop("'details' must be type logical.", call. = FALSE)
  }

  #--- Extract values from mraster ---#

  vals <- terra::values(mraster)

  #--- Determine index of each cell so to map values correctly without NA ---#

  idx <- which(complete.cases(vals))

  valsOut <- vals

  #--- conduct unsupervised k-means with center/scale parameters based on algorithm ---#

  message("K-means being performed on ", terra::nlyr(mraster), " layers with ", nStrata, " centers.")

  km_clust <- stats::kmeans(scale(vals[idx], center = center, scale = scale), centers = nStrata, iter.max = iter, algorithm = algorithm)

  #--- convert k-means values back to original mraster extent ---#
  valsOut[idx] <- km_clust$cluster

  kmv <- suppressWarnings(terra::setValues(mraster[[1]], valsOut))
  names(kmv) <- "strata"

  #--- plot if requested ---#

  if (isTRUE(plot)) {
    terra::plot(kmv, main = "K-means clusters", type = "classes")
  }

  #--- write file to disc ---#

  if (!is.null(filename)) {
    terra::writeRaster(x = kmv, filename = filename, overwrite = overwrite)
    message("Output raster written to disc.")
  }

  #--- Output based on 'details' to return raster alone or list with details ---#

  if (isTRUE(details)) {

    #--- create list to assign k-means info and output raster ---#

    out <- list(details = km_clust, raster = kmv)

    return(out)
  } else {

    #--- just output raster ---#

    return(kmv)
  }
}
