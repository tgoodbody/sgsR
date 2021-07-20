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
#' @param scale Logical. Value indicating whether the variables should be scaled to have unit variance
#' @param plot Logical. Plots output strata raster and visualized
#'  strata with boundary dividers.
#' @param details Logical. If \code{FALSE} (default) output is only
#'  stratification raster. If \code{TRUE} return a list
#' where \code{$details} is additional stratification information and
#'  \code{$raster} is the output stratification spatRaster.
#'  @param ... Additional arguments to be passed to \code{\link[stats]{kmeans}} function.
#'
#' @return output stratification \code{spatRaster}, or a list when \code{details = TRUE}.
#' 
#' @examples 
#' #--- Load raster and access files ---#
#' r <- system.file("extdata","wall_metrics_small.tif", package = "sgsR"
#' mr <- terra::rast(r)
#' 
#' #--- perform stratification using k-means ---#
#' kmeans <- strat_kmeans(mraster = mr, nStrata = 5)
#' 
#' kmeans <- strat_kmeans(mraster = mr, nStrata = 5, iter = 1000, plot = TRUE, details = TRUE)
#' 
#' kmeans <- strat_kmeans(mraster = mr, nStrata = 5, iter = 1000, plot = TRUE, tempfile(fileext = ".tif"), overwrite = TRUE)
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
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  }

  if (!is.numeric(nStrata)) {
    stop("'nStrata' must be type numeric")
  }

  if (!is.numeric(iter)) {
    stop("'iter.max' must be type numeric")
  }

  if (!is.character(algorithm)) {
    stop("'algorithm' must be type character")
  }

  if (algorithm != "Lloyd" && algorithm != "MacQueen") {
    stop("Unknown algorithm '", algorithm, "' selected. Please use 'Lloyd' (default) or 'MacQueen'")
  }

  if (!is.logical(center)) {
    stop("'center' must be type logical")
  }

  if (!is.logical(scale)) {
    stop("'scale' must be type logical")
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical")
  }

  if (!is.logical(details)) {
    stop("'details' must be type logical")
  }


  #--- Extract values from mraster ---#

  vals <- terra::values(mraster)

  #--- Determine index of each cell so to map values correctly without NA ---#

  vals[!is.finite(vals)] <- NA

  #--- conduct unsupervised k-means with center/scale parameters based on algorithm ---#

  message("K-means being performed on ", terra::nlyr(mraster), " layers with ", nStrata, " centers.")

  km_clust <- stats::kmeans(scale(na.omit(vals), center = center, scale = scale), centers = nStrata, iter.max = iter, algorithm = algorithm, ...)

  #--- convert k-means values back to original mraster extent ---#

  vals[is.finite(vals)] <- km_clust$cluster

  kmv <- terra::setValues(mraster[[1]], vals)
  names(kmv) <- "strata"


  #--- plot if requested ---#

  if (isTRUE(plot)) {
    terra::plot(kmv, main = "K-means clusters", type = "classes")
  }

  #--- write file to disc ---#

  if (!is.null(filename)) {
    terra::writeRaster(kmv, filename, overwrite = overwrite)
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
