#' Raster principal components
#'
#' @description Calculate and rasterize principal components from a metric raster
#'
#' @family calculate functions
#'
#' @inheritParams sample_srs
#' @inheritParams strat_kmeans
#'
#' @param nComp Numeric. Value indicating number of principal components to be rasterized.
#' prior to analysis.
#' @param ... Additional arguments to be passed to \code{\link[stats]{prcomp}}.
#'
#' @importFrom stats na.exclude na.omit prcomp
#'
#' @return Output raster with specified number of principal components as layers.
#'
#' @examples
#' #--- Load raster ---#
#' r <- system.file("extdata", "wall_metrics.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' calculate_pcomp(
#'   mraster = mr,
#'   nComp = 2,
#'   plot = TRUE
#' )
#'
#' pcomp <- calculate_pcomp(
#'   mraster = mr,
#'   nComp = 3,
#'   details = TRUE
#' )
#'
#' #--- Display principal component details ---#
#' pcomp$pca
#'
#' #--- Display importance of components ---#
#' summary(pcomp$pca)
#' @author Tristan R.H. Goodbody
#'
#' @export

calculate_pcomp <- function(mraster,
                            nComp,
                            center = TRUE,
                            scale = TRUE,
                            plot = FALSE,
                            details = FALSE,
                            filename = NULL,
                            overwrite = FALSE,
                            ...) {

  #--- set global vars ---#

  rout <- NULL

  #--- error handling ---#

  if (!inherits(mraster, "SpatRaster")) {
    stop("'raster' must be type SpatRaster", call. = FALSE)
  }

  if (!is.numeric(nComp)) {
    stop("'ncomp' must be type numeric")
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

  if (nComp > terra::nlyr(mraster)) {
    stop("nComp must be <= the number of layers in 'mraster'")
  }

  #--- Extract values from mraster ---#

  vals <- terra::values(mraster)

  #--- perform PCA ---#

  PCA <- stats::prcomp(
    formula = ~.,
    data = as.data.frame(vals),
    center = center,
    scale. = scale,
    na.action = na.exclude,
    ...
  )

  #--- extract cell level pca values ---#
  pcavals <- as.data.frame(PCA$x)

  #--- apply cell-level allocation of values---#

  rs <- apply(X = pcavals, MARGIN = 2, FUN = terra::setValues, x = mraster[[1]])

  #--- stack pca rasters ---#

  pcaRout <- terra::rast(rs[1:nComp])

  #--- rename ---#

  names(pcaRout) <- rep(paste0("PC", seq(1, nComp, 1)))

  #--- write file to disc ---#

  if (!is.null(filename)) {
    terra::writeRaster(rout, filename, overwrite = overwrite)
  }

  if (isTRUE(plot)) {

    #--- Plot components ---#
    terra::plot(pcaRout)
  }

  if (isTRUE(details)) {

    #--- create list to assign pca info and output raster ---#

    out <- list(pca = PCA, raster = pcaRout)

    return(out)
  } else {

    #--- just output raster ---#

    return(pcaRout)
  }
}
