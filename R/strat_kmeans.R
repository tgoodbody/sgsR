#' k-means stratification
#'
#' @description Stratify metrics raster using k-means algorithm
#' @family stratify functions
#'
#' @inheritParams strat_breaks
#'
#' @param mraster spatRaster. ALS metrics raster.
#' @param nstrata Character. Number of desired strata.
#' @param iter.max Numeric. The maximum number of iterations allowed.
#' @param algorithm Character. \code{Lloyd} (default) or
#' \code{MacQueen} algorithms.
#' @param scale Logical. Determines whether scaling and centering
#'  of data should be conducted prior to analysis.
#' @param plot Logical. Plots output strata raster and visualized
#'  strata with boundary dividers.
#' @param details Logical. If \code{FALSE} (default) output is only
#'  stratification raster. If \code{TRUE} return a list
#' where \code{$details} is additional stratification information and
#'  \code{$raster} is the output stratification spatRaster.
#'
#' @importFrom magrittr %>%
#' @importFrom methods is
#'
#' @return output stratification \code{spatRaster}, or a list when \code{details = TRUE}.
#'
#' @export

strat_kmeans <- function(mraster,
                         nstrata,
                         iter.max = 500,
                         algorithm = "Lloyd",
                         scale = TRUE,
                         plot = FALSE,
                         details = FALSE,
                         filename = NULL,
                         ...) {

  #--- Error management ---#

  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  }

  if (!is.numeric(nstrata)) {
    stop("'nstrata' must be type numeric")
  }

  if (!is.numeric(iter.max)) {
    stop("'iter.max' must be type numeric")
  }

  if (!is.character(algorithm)) {
    stop("'algorithm' must be type character")
  }

  if (algorithm != "Lloyd" && algorithm != "MacQueen") {
    stop("Unknown algorithm '", algorithm, "' selected. Please use 'Lloyd' (default) or 'MacQueen'")
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

  idx <- !is.na(vals)

  #--- conduct unsupervised k-means with center/scale parameters based on algorithm ---#

  if (algorithm == "Lloyd") {
    if (isTRUE(scale)) {
      message("K-means being performed on ", terra::nlyr(mraster), " layers with ", nstrata, " centers. Data have been centered and scaled")

      km_clust <- stats::kmeans(scale(vals[idx], center = TRUE, scale = TRUE), centers = nstrata, iter.max = iter.max, algorithm = algorithm)
    } else {
      message("K-means being performed on ", terra::nlyr(mraster), " layers with ", nstrata, " centers, Data have not been centered or scaled.")

      km_clust <- stats::kmeans(scale(vals[idx], center = FALSE, scale = FALSE), centers = nstrata, iter.max = iter.max, algorithm = algorithm)
    }
  }

  #--- conduct unsupervised k-means with center/scale parameters based on algorithm ---#

  else if (algorithm == "MacQueen") {
    if (isTRUE(scale)) {
      message("K-means being performed on ", terra::nlyr(mraster), " layers with ", nstrata, " centers. Data have been centered and scaled")

      km_clust <- stats::kmeans(scale(vals[idx], center = TRUE, scale = TRUE), centers = nstrata, iter.max = iter.max, algorithm = algorithm)
    } else {
      message("K-means being performed on ", terra::nlyr(mraster), " layers with ", nstrata, " centers, Data have not been centered or scaled.")

      km_clust <- stats::kmeans(scale(vals[idx], center = FALSE, scale = FALSE), centers = nstrata, iter.max = iter.max, algorithm = algorithm)
    }
  }

  #--- convert k-means values back to original mraster extent ---#

  vals[idx] <- km_clust$cluster

  kmv <- terra::setValues(mraster[[1]], vals)
  names(kmv) <- "strata"



  #--- plot if requested ---#

  if (isTRUE(plot)) {

    #--- make plot using diverging colour palette ---#

    ncols <- nstrata
    col <- RColorBrewer::brewer.pal(ncols, "Set3")

    terra::plot(kmv, main = "K-means clusters", col = col, type = "classes")
  }

  #--- write file to disc ---#

  if (!is.null(filename)) {
    writeRaster(kmv, filename, overwrite = TRUE, ...)
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
