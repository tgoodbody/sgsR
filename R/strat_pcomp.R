#' Principal components stratification
#'
#' @description Stratify metrics raster using principal components and quantile breaks
#' @family stratify functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams strat_breaks
#' @inheritParams calculate_pcomp
#' @param mraster Spatraster. Covariate raster to stratify.
#' @param nStrata Numeric. Number of strata for primary principal component.
#' @param nStrata2 Numeric. Number of strata for secondary principal component.
#' @param ... Additional arguments to be passed to \code{\link[stats]{prcomp}}.
#'
#' @return Returns an output stratification \code{spatRaster} or a list when \code{details = TRUE}.
#'
#' When a list is returned:
#' \enumerate{
#' \item \code{details} is a list output of the \code{\link[stats]{prcomp}} function
#' \item \code{raster} is a stratified \code{spatRaster} based on the PCA
#' \item \code{plot} is a \code{ggplot} scatter plot object where strata are colour coded
#' and strata boundaries are delineated
#' }
#'
#' @examples
#' #--- Load raster and access files ---#
#' r <- system.file("extdata", "wall_metrics.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' strat_pcomp(
#'   mraster = mr,
#'   nStrata = 5,
#'   plot = TRUE
#' )
#'
#' strat_pcomp(
#'   mraster = mr,
#'   nStrata = 4,
#'   nStrata2 = 4,
#'   plot = TRUE,
#'   details = TRUE
#' )
#'
#' strat_pcomp(
#'   mraster = mr,
#'   nStrata = 3,
#'   nStrata2 = 3,
#'   filename = tempfile(fileext = ".tif")
#' )
#' @author Tristan R.H. Goodbody
#'
#' @export

strat_pcomp <- function(mraster,
                        nStrata,
                        nStrata2 = NULL,
                        center = TRUE,
                        scale = TRUE,
                        plot = FALSE,
                        samp = 1,
                        details = FALSE,
                        filename = NULL,
                        overwrite = FALSE,
                        ...) {

  #--- set global vars ---#

  raster <- class1 <- class2 <- PC1 <- PC2 <- NULL

  #--- Error management ---#

  if (!inherits(mraster, "SpatRaster")) {
    stop("mraster must be type SpatRaster", call. = FALSE)
  }

  if (!is.numeric(nStrata)) {
    stop("'nStrata' must be type numeric")
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

  if (!is.numeric(samp)) {
    stop("'samp' must be type numeric")
  }

  if (!is.logical(details)) {
    stop("'details' must be type logical")
  }

  #--- Extract values from mraster ---#

  vals <- terra::values(mraster)


  if (is.null(nStrata2)) {

    #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#

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

    pcagrps <- pcavals %>%
      #--- define nStrata classes ---#
      dplyr::mutate(class = dplyr::ntile(PC1, nStrata))


    #--- set newly stratified values ---#

    rout <- terra::setValues(mraster[[1]], pcagrps$class)
    names(rout) <- "strata"
  } else {
    if (!is.numeric(nStrata2)) {
      stop("'nStrata2' must be type numeric")
    }

    #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#

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

    #--- Split PCA distribution in to number specified by 'nStrata' ---#

    pcagrps <- pcavals %>%
      #--- define nStrata classes ---#
      dplyr::mutate(class1 = dplyr::ntile(PC1, nStrata)) %>%
      #--- group by class to sub stratify ---#
      dplyr::group_by(class1) %>%
      #--- define nStrata2 classes ---#
      dplyr::mutate(class2 = dplyr::ntile(PC2, nStrata2)) %>%
      #--- combine classes ---#
      dplyr::group_by(class1, class2) %>%
      #--- establish newly formed unique class ---#
      dplyr::mutate(class = dplyr::cur_group_id()) %>%
      #--- ensure NA's are transfered ---#
      dplyr::mutate(class = ifelse(is.na(class1), NA, class))


    #--- set newly stratified values ---#

    rout <- terra::setValues(mraster[[1]], pcagrps$class)
    names(rout) <- "strata"
  }

  if (isTRUE(plot)) {
    if (samp > 1 | samp < 0) {
      stop("'samp' must be > 0 & <= 1")
    }

    terra::plot(rout)
  }

  #--- write file to disc ---#

  if (!is.null(filename)) {
    terra::writeRaster(rout, filename, overwrite = overwrite)
  }

  #--- Output based on 'details' to return raster alone or list with details ---#

  if (isTRUE(details)) {

    #--- create classplot summary ---#

    coordsgrps <- pcagrps %>%
      dplyr::group_by(class) %>%
      dplyr::arrange(class) %>%
      stats::na.omit() %>%
      tidyr::nest() %>%
      dplyr::ungroup()

    p <- classPlot(as.data.frame(pcagrps),
      coordsgrps,
      mraster = PC1,
      mraster2 = PC2,
      samp
    )

    #--- create list to assign pca info and output raster ---#

    out <- list(details = PCA, raster = rout, plot = p)

    return(out)
  } else {

    #--- just output raster ---#

    return(rout)
  }
}
