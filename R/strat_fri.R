#' Stratify Forest Resources Inventory (FRI)
#'
#' @description Stratify an input forest resources inventory (FRI) polygon.
#'
#' @family stratify functions
#'
#' @inheritParams strat_breaks
#' @param fri sf. Forest resources inventory polygon coverage.
#' @param attribute Character. Name of attribute within \code{fri} that will be stratified
#' @param features vector/list. Vector or list of vectors of features within \code{attribute}
#' to guide stratification.
#' @param raster spatRaster. Raster for polygon to raster conversion. If \code{map = TRUE}
#' raster must contain a layer named "strata".
#' @param plot Logical. Plots output spatRaster.
#' @param details Logical. If \code{FALSE} (default) output is spatRaster object of
#' stratified forest resources inventory attributes. If \code{TRUE} return a list
#' where \code{$outRaster} is the stratified forest resources inventory attributes,
#' \code{$lookUp} is the lookup table for the stratification, and \code{fripoly} is the
#' forest resources inventory poly with \code{attribute} and corresponding \code{strata}
#'
#' @examples
#' #--- load input metrics raster ---#
#' raster <- system.file("extdata","kmeans.tif", package = "sgsR")
#' srasterkmeans <- terra::rast(raster)
#' 
#' #--- read polygon coverage ---#
#' poly <- system.file("extdata","inventory_polygons.shp", package = "sgsR")
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
#' features <- c("poor","rich","medium")
#' 
#' srasterfri <- strat_fri(fri = fri, 
#'                         attribute = attribute, 
#'                         features = features, 
#'                         raster = mraster, 
#'                         plot = TRUE)
#'                         
#' #--- or as multiple lists ---#
#' 
#' g1 <- "poor"
#' g2 <- c("rich", "medium")
#' 
#' features <- list(g1,g2)
#' 
#' srasterfri <- strat_fri(fri = fri, 
#'                         attribute = attribute, 
#'                         features = features, 
#'                         raster = mraster,
#'                         stack = TRUE, 
#'                         plot = TRUE,
#'                         details = TRUE)
#' 
#'
#' @author Tristan R.H. Goodbody
#'
#' @importFrom methods is
#'
#' @return A spatRaster object.
#' @export

strat_fri <- function(fri,
                      attribute,
                      features,
                      raster,
                      filename = NULL,
                      overwrite = FALSE,
                      plot = FALSE,
                      details = FALSE,
                      ...
){
  
  #--- error handling ---#
  
  if (!inherits(fri, "sf")) {
    stop("'access' must be an 'sf' object")
  }
  
  if (!inherits(sf::st_geometry(fri), "sfc_POLYGON") && !inherits(sf::st_geometry(fri), "sfc_MULTIPOLYGON")) {
    stop("'fri' geometry type must be 'POLYGON' or 'MULTIPOLYGON'")
  }
  
  if (!is.character(attribute)) {
    stop("'attribute' must be type character")
  }
  
  if(!is.vector(features) && !is.list(features)){
    stop("'features' must supply a vector or a list of vectors")
  }
  
  if (!inherits(raster, "SpatRaster")) {
    stop("'raster' must be type SpatRaster", call. = FALSE)
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
  
  if (any(!glue::glue('{attribute}') %in% names(fri))) {
    stop(glue('fri does not have a layer named {attribute}'))
  }
  
  #--- check that features are not duplicated across proposed attribute classes ---#
  
  unFeat <- unlist(features)
  
  unFeatObjs <- unFeat[duplicated(unFeat)]
  
  if(length(unFeatObjs) > 0){
    stop(glue::glue("{unFeatObjs*}", .transformer = collapse_transformer(sep = ", ", last = " and ")), ". Are duplicated in 'features'.")
    
  }
  
  #--- begin fri polygon manipulation ---#
  
  fri <- fri %>% 
    dplyr::select(glue::glue('{attribute}'))
  
  #--- generate number of groups ---#

  strata <- as.list(seq(1,length(features),1))
  
  #--- create lookup table for values and associated strataes ---#
  
  lookUp <- do.call(rbind, Map(data.frame, strata=strata, features=features))
  
  #--- check if features specifiec are included within fri attribute ---#
  
  if (any(!lookUp$features %in% fri[[1]])) {
    stop("'attribute' does not have specified 'features'.")
  }

  #--- create new column with mutated values based on features and associated group ---#
  
  fripoly <- fri %>% 
    dplyr::mutate(strata = dplyr::case_when(!!!rlang::parse_exprs(glue::glue('{attribute} %in% "{lookUp$features}" ~ "{lookUp$strata}"')))) %>%
    na.omit() %>%
    terra::vect()
  
  #--- define raster ---#
  raster <- raster[[1]]
  
  #--- rasterize vector ---#
  outfri <- terra::rasterize(x = fripoly, y = raster, field = "strata")
  
  #--- write file to disc ---#
  
  if (!is.null(filename)) {
    terra::writeRaster(outfri, filename, overwrite = overwrite, ...)
  }
  
  #--- plot if requested
  
  if(isTRUE(plot)){
    
    terra::plot(outfri)
    
  }
  
  #--- output details if desired ---#
  
  if (isTRUE(details)) {
    
    #--- output metrics details along with stratification raster ---#
    
    output <- list(outRaster = outfri, lookUp = lookUp, fripoly = fripoly)
    
    #--- output samples dataframe ---#
    
    return(output)
    
  } else {
    
    #--- just output raster ---#
    
    return(outfri)
  }
}

#--- glue transformer ---#

collapse_transformer <- function(regex = "[*]$", ...) {
  function(text, envir) {
    collapse <- grepl(regex, text)
    if (collapse) {
      text <- sub(regex, "", text)
    }
    res <- identity_transformer(text, envir)
    if (collapse) {
      glue_collapse(res, ...)  
    } else {
      res
    }
  }
}