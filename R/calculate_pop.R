#' Population descriptors
#'
#' @description Population matrices and descriptions of metric raster data
#'
#' @family calculate functions
#'
#' @description Calculates population level statistics including principal components, quantile matrix, and covariance matrix
#' needed necessary for \code{\link{calculate_lhsOpt}}. Outputs can also be used as an input for \code{\link{sample_ahels}}.
#'
#' @inheritParams strat_kmeans
#'
#' @param PCA Logical. Calculates principal component loadings of the population for PCA similarity factor testing.
#' \code{default = FALSE}.
#' @param matQ Logical. Calculates quantile matrix of the population for quantile comparison testing.
#' \code{default = TRUE}.
#' @param nQuant Numeric. Number of quantiles to divide the population into for \code{matQ}.
#' \code{default = 10}.
#' @param matCov Logical. Calculates covariate matrix of the population. Needed for Kullback–Leibler divergence testing.
#' \code{default = TRUE}. Requires \code{matQ = TRUE}.
#'
#' @references
#' Malone BP, Minansy B, Brungard C. 2019. Some methods to improve the utility of conditioned Latin hypercube sampling. PeerJ 7:e6451 DOI 10.7717/peerj.6451
#'
#' @return List of matrices to be used as input for \code{\link{calculate_lhsOpt}}.
#'
#' @examples
#' #--- Load raster and access files ---#
#' r <- system.file("extdata", "mraster.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' calculate_pop(mraster = mr)
#'
#' @note
#' Special thanks to Dr. Brendan Malone for the original implementation of this algorithm.
#'
#' @author Tristan R.H. Goodbody
#'
#' @export

calculate_pop <- function(mraster,
                          PCA = FALSE,
                          matQ = TRUE,
                          nQuant = 10,
                          matCov = TRUE) {
  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster.", call. = FALSE)
  }

  if (!is.logical(PCA)) {
    stop("'PCA' must be type logical.", call. = FALSE)
  }

  if (!is.logical(matQ)) {
    stop("'matQ' must be type logical.", call. = FALSE)
  }

  if (!is.logical(matCov)) {
    stop("'matCov' must be type logical.", call. = FALSE)
  }

  #--- determine number of bands in 'mraster' ---#

  nb <- terra::nlyr(mraster)

  #--- Extract values from mraster ---#

  vals <- terra::values(mraster)

  #--- Determine index of each cell so to map values correctly without NA ---#

  vals[!is.finite(vals)] <- NA

  #--- Remove NA / NaN / Inf values ---#

  vals <- vals %>%
    as.data.frame() %>%
    dplyr::filter(stats::complete.cases(.))

  #--- PCA loadings for the population ---#

  if (isTRUE(PCA)) {
    #--- perform PCA analysis for the population to determine variance in each component ---#

    pca <- stats::prcomp(vals, scale = TRUE, center = TRUE)

    #--- extract pca scores ---#

    pcaScores <- as.data.frame(pca$x)

    #--- extract PCA loadings ---#

    pcaLoad <- matrix(NA, ncol = nb, nrow = nb)

    for (i in 1:nb) {
      pcaLoad[i, ] <- as.matrix(t(pca$rotation[i, ]))
    }
  } else {
    pcaLoad <- NULL
  }

  #--- Quantiles of the population ---#

  if (isTRUE(matQ)) {
    if (!is.numeric(nQuant)) {
      stop("'nQuant' must be type numeric.", call. = FALSE)
    }

    matQ <- mat_quant(
      vals,
      nQuant,
      nb
    )
  } else {
    matQ <- NULL
  }

  #--- covariate hypercube for KL divergence test ---#

  if (isTRUE(matCov)) {
    if (is.null(matQ)) {
      stop("Covariance matrix creation requires quantile matrix. Set 'matQ = TRUE'.", call. = FALSE)
    }

    #--- create covariate matrix of the quantiles ---#

    message("Creating covariance matrix.")

    matCov <- mat_cov(vals, nQuant, nb, matQ)
  } else {
    matCov <- NULL
  }

  #--- create list of outputs ---#

  lout <- list(
    values = vals,
    pcaLoad = pcaLoad,
    matQ = matQ,
    matCov = matCov
  )

  #--- remove NULL list objects ---#

  lout <- lout[!sapply(lout, is.null)]

  #--- output list ---#

  return(lout)
}
