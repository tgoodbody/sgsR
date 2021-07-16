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

#' @importFrom methods is
#'
#' @return Output raster with specified number of principal components as layers.
#'
#' @export

calculate_pcomp <- function(mraster = NULL,
                            nComp = NULL,
                            center = TRUE,
                            scale = TRUE,
                            plot = FALSE,
                            details = FALSE,
                            ...) {

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

  #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#

  PCA <- prcomp(
    formula = ~.,
    data = as.data.frame(vals),
    center = center,
    scale. = scale,
    na.action = na.exclude,
    ...
  )

  #--- extract cell level pca values ---#
  pcavals <- as.data.frame(PCA$x)

  #--- create loop to allocate values to cell level taking into account potential NA using 'idx' ---#
  rs <- list()

  for (i in 1:nComp) {
    rs[[i]] <- terra::setValues(mraster[[1]], pcavals[, i, drop = FALSE])
  }

  #--- stack pca rasters ---#

  pcaRout <- terra::rast(rs)

  #--- rename ---#

  names(pcaRout) <- rep(paste0("PC", seq(1, nComp, 1)))

  if (isTRUE(plot)) {

    #--- visualize scree plot ---#
    terra::plot(pcaRout[[1:nComp]])
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
