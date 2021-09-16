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
#' @param map Logical. Default = \code{FALSE}. If \code{TRUE}, map values from fristrata 
#' with strata from \code{raster} to create secondary strata tier.
#' @param stack Logical. Default = \code{FALSE}. If \code{TRUE}, output raster will be
#' 3 layers: \code{fristrata, rasterstrata, frirasterStrata}.
#' @param plot Logical. Plots output spatRaster.
#' @param details Logical. If \code{FALSE} (default) output is spatRaster object of
#' stratified forest resources inventory attributes. If \code{TRUE} return a list
#' where \code{$outRaster} is the stratified forest resources inventory attributes,
#' \code{$lookUp} is the lookup table for the stratification, and \code{fripoly} is the
#' forest resources inventory poly with \code{attribute} and corresponding \code{class}
#'
#' @importFrom methods is
#'
#' @return A spatRaster object.


strat_fri <- function(fri,
                      attribute,
                      features,
                      raster,
                      map = FALSE,
                      stack = FALSE,
                      filename = NULL,
                      overwrite = FALSE,
                      plot = FALSE,
                      details = FALSE
){
  
  #--- error handling ---#
  
  if (!inherits(fri, "sf")) {
    stop("'access' must be an 'sf' object")
  }
  
  if (!inherits(sf::st_geometry(fri), "sfc_POLYGON") && !inherits(sf::st_geometry(fri), "sfc_MULTIPOLYGON")) {
    stop("'fri' geometry type must be 'POLYGON' or 'MULTIPOLYGON'")
  }
  
  if (!inherits(raster, "SpatRaster")) {
    stop("'raster' must be type SpatRaster", call. = FALSE)
  }
  
  if (!is.character(attribute)) {
    stop("'attribute' must be type character")
  }
  
  if(!is.vector(features) && !is.list(features)){
    stop("'features' must supply a vector or a list of vectors")
  }
  
  if (!is.logical(map)) {
    stop("'map' must be type logical")
  }
  
  if (!is.logical(stack)) {
    stop("'stack' must be type logical")
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
  
  fri <- fri %>% 
    dplyr::select(glue::glue('{attribute}'))
  
  #--- generate number of groups ---#

  class <- as.list(seq(1,length(features),1))
  
  #--- create lookup table for values and associated classes ---#
  
  lookUp <- do.call(rbind, Map(data.frame, class=class, features=features))
  
  #--- check if features specifiec are included within fri attribute ---#
  
  if (any(!lookUp$features %in% fri[[1]])) {
    stop("'attribute' does not have specified 'features'.")
  }

  #--- create new column with mutated values based on features and associated group ---#
  
  fripoly <- fri %>% 
    dplyr::mutate(class = dplyr::case_when(!!!rlang::parse_exprs(glue::glue('{attribute} %in% "{lookUp$features}" ~ "{lookUp$class}"')))) %>%
    na.omit() %>%
    terra::vect()
  
  #--- define raster ---#
  raster <- raster[[1]]
  
  #--- rasterize vector ---#
  outfri <- terra::rasterize(x = fripoly, y = raster, field = "class")
  
  #--- combine fri raster with stratification raster ---#
  
  if(isTRUE(map)){
    
    message("fri stratification complete. Mapping stratified fri with 'raster' strata")
    
    if (any(!c("strata") %in% names(raster))) {
      stop("'map = TRUE' but 'raster' does not have a layer named 'strata'")
    }
    
    joined <- c(outfri,raster[["strata"]])
    
    featuresJoin <- terra::values(joined)
    
    oclass <- featuresJoin %>%
      as.data.frame() %>%
      dplyr::group_by(class, strata) %>%
      #--- establish newly formed unique class ---#
      dplyr::mutate(strata_combined = paste0(class,strata)) %>%
      #--- ensure NA's are transfered ---#
      dplyr::mutate(strata_combined = ifelse(is.na(class) | is.na(strata), NA, strata_combined)) %>%
      dplyr::rename(fristrata = class,
                    rasterstrata = strata,
                    frirasterStrata = strata_combined)
    
    #--- create lookUp table ---#
    
    lookUp <- distinct(oclass) %>% 
      na.omit() %>%
      as.data.frame()
    
    #--- set newly stratified values ---#
    
    outfrijoin <- terra::setValues(raster, oclass$frirasterStrata)
    names(outfri) <- "frirasterStrata"
    
  }
  
  if(isTRUE(stack)){
    
    message("Stacking fri strata, raster strata, and their combination.")
    
    #--- stack 3 rasters if requested ---#
    
    outfri <- c(outfri,raster,outfrijoin)
    names(outfri) <- c("fristrata","rasterstrata","frirasterStrata")
    
  }
  
  #--- if not stacking rename for output ---#
  
  if(exists("outfrijoin")){

    outfri <- outfrijoin
    
  }
  
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