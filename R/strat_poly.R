#' Stratify using polygons
#'
#' @description Stratify based on polygon coverage attributes and features.
#'
#' @family stratify functions
#'
#' @inheritParams strat_breaks
#' @param poly sf. Input polygon coverage. e.g. - forest resources inventory coverage.
#' @param attribute Character. Name of attribute within \code{poly} that will be stratified
#' @param features vector/list. Vector or list of vectors of features within \code{attribute}
#' to guide stratification.
#' @param raster spatRaster. Raster template to enable polygon to raster conversion.
#' @param plot Logical. Plots output spatRaster.
#' @param details Logical. If \code{FALSE} (default) output is spatRaster object of
#' stratified polygon attributes. If \code{TRUE} return a list
#' where \code{$outRaster} is the stratified attributes, \code{$lookUp} is the lookup table
#' for the stratification, and \code{poly} is the defined polygon \code{attribute} with corresponding
#' \code{features / strata}
#'
#' @examples
#' #--- load input metrics raster ---#
#' raster <- system.file("extdata", "kmeans.tif", package = "sgsR")
#' sraster <- terra::rast(raster)
#'
#' #--- read polygon coverage ---#
#' poly <- system.file("extdata", "inventory_polygons.shp", package = "sgsR")
#' fri <- sf::st_read(poly)
#'
#' #--- stratify polygon coverage ---#
#' #--- specify polygon attribute to stratify ---#
#'
#' attribute <- "NUTRIENTS"
#'
#' #--- specify features within attribute & how they should be grouped ---#
#' #--- as a single vector ---#
#'
#' features <- c("poor", "rich", "medium")
#'
#' srasterpoly <- strat_poly(
#'   poly = fri,
#'   attribute = attribute,
#'   features = features,
#'   raster = sraster,
#'   plot = TRUE
#' )
#'
#' #--- or as multiple lists ---#
#'
#' g1 <- "poor"
#' g2 <- c("rich", "medium")
#'
#' features <- list(g1, g2)
#'
#' srasterpoly <- strat_poly(
#'   poly = fri,
#'   attribute = attribute,
#'   features = features,
#'   raster = sraster,
#'   plot = TRUE,
#'   details = TRUE
#' )
#' @author Tristan R.H. Goodbody
#'
#' @importFrom methods is
#'
#' @return A spatRaster object.
#' @export

strat_poly <- function(poly,
                       attribute,
                       features,
                       raster,
                       filename = NULL,
                       overwrite = FALSE,
                       plot = FALSE,
                       details = FALSE,
                       ...) {

  #--- error handling ---#

  if (!inherits(poly, "sf")) {
    stop("'poly' must be an 'sf' object")
  }

  if (!inherits(sf::st_geometry(poly), "sfc_POLYGON") && !inherits(sf::st_geometry(poly), "sfc_MULTIPOLYGON")) {
    stop("'poly' geometry type must be 'POLYGON' or 'MULTIPOLYGON'")
  }

  if (!is.character(attribute)) {
    stop("'attribute' must be type character")
  }

  if (!is.vector(features) && !is.list(features)) {
    stop("'features' must supply a vector or a list of vectors")
  }

  if (!inherits(raster, "SpatRaster")) {
    stop("'raster' must be type SpatRaster")
  }

  if (!is.logical(overwrite)) {
    stop("'overwrite' must be either TRUE or FALSE")
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical")
  }

  if (!is.logical(details)) {
    stop("'details' must be type logical")
  }

  #--- subset inventory polygon ---#

  if (!any(grepl(attribute, names(poly)))) {
    stop(paste0("'poly' does not have a layer named ", attribute))
  }

  #--- check that features are not duplicated across proposed attribute classes ---#

  unFeat <- unlist(features)

  if (any(is.na(unFeat))) {
    message("'features' contains NA. Is this on purpose?")
  }

  unFeatObjs <- unFeat[duplicated(unFeat)]

  if (length(unFeatObjs) > 0) {
    stop(paste(c("Repeated within 'features':", unFeatObjs), collapse = " "))
  }

  #--- begin polygon manipulation ---#

  poly <- poly %>%
    dplyr::select(!!attribute)

  #--- generate number of groups ---#

  strata <- as.list(seq(1, length(features), 1))

  #--- create lookup table for values and associated strataes ---#

  lookUp <- do.call(rbind, Map(data.frame, strata = strata, features = features))

  #--- check if features specifiec are included within poly attribute ---#

  if (any(!lookUp$features %in% poly[[1]])) {
    stop("'attribute' does not have specified 'features'.")
  }

  #--- create new column with mutated values based on features and associated group ---#

  poly <- poly %>%
    dplyr::rename(features = attribute) %>%
    dplyr::left_join(lookUp, by = "features") %>%
    na.omit() %>%
    terra::vect()

  #--- ensure that strata values are integer for mapping ---#

  poly$strata <- as.integer(poly$strata)

  #--- rasterize vector ---#
  outpolyrast <- terra::rasterize(x = poly, y = raster[[1]], field = "strata")

  #--- write file to disc if requested ---#

  if (!is.null(filename)) {
    terra::writeRaster(x = outpolyrast, filename = filename, overwrite = overwrite, ...)
  }

  #--- plot if requested ---#

  if (isTRUE(plot)) {
    terra::plot(outpolyrast)
  }

  #--- output details if requested ---#

  if (isTRUE(details)) {

    #--- output metrics details along with stratification raster ---#

    output <- list(outRaster = outpolyrast, lookUp = lookUp, poly = poly)

    #--- output samples dataframe ---#

    return(output)
  } else {

    #--- just output raster ---#

    return(outpolyrast)
  }
}