#' Map 2 stratified rasters
#'
#' @description Map stratified rasters to a combined stratification.
#'
#' @family stratify functions
#'
#' @inheritParams strat_breaks
#' @inheritParams strat_fri
#' @param sraster spatRaster. Primary stratification raster.
#' @param sraster2 spatRaster. Secondary stratification raster.
#' @param stack Logical. Default = \code{FALSE}. If \code{TRUE}, output raster will be
#' 3 layers: \code{strata, strata2, stratamapped}.
#' @param details Logical. If \code{FALSE} (default) output is a mapped stratified spatRaster object.
#' If \code{TRUE} return a list where \code{$outRaster} is the mapped stratified raster, and
#' \code{$lookUp} is the lookup table for the stratification.
#' 
#' @section Mapping:
#' The mapping algorithm will take the stratification from \code{sraster} and combine it with
#' overlying strata values in \code{sraster2}. This will result in a \code{stratamapped} attribute
#' where the values from both inputs are combined. 
#' 
#' i.e.
#' 
#' If \code{strata = 1} and \code{strata2 = 1} then \code{stratamapped = 11}.
#' 
#' If \code{strata = 2} and \code{strata2 = 14} then \code{stratamapped = 214}.
#'
#'@examples
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
#' #--- map srasters ---#
#' strat_map(sraster = srasterfri,
#'           sraster2 = srasterkmeans,
#'           plot = TRUE)
#' 
#' strat_map(sraster = srasterfri,
#'           sraster2 = srasterkmeans,
#'           stack = TRUE,
#'           details = TRUE,
#'           plot = TRUE)
#'           
#' @importFrom methods is
#'
#' @return A spatRaster object.
#' @export


strat_map <- function(sraster,
                      sraster2,
                      stack = FALSE,
                      filename = NULL,
                      overwrite = FALSE,
                      plot = FALSE,
                      details = FALSE,
                      ...){
  
  #--- error handling ---#
  
  if (!inherits(sraster, "SpatRaster")) {
    stop("'sraster' must be type SpatRaster", call. = FALSE)
  }
  
  if (!inherits(sraster2, "SpatRaster")) {
    stop("'sraster' must be type SpatRaster", call. = FALSE)
  }
  
  if (!is.logical(stack)) {
    stop("'stack' must be type logical")
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
  
  #--- error handling for raster inputs ---#
  
  if(terra::nlyr(sraster) > 1){
    stop("sraster mustonly contain 1 layer. Please subset the layer you would like to use for mapping.")
  }
  
  if(terra::nlyr(sraster2) > 1){
    stop("sraster2 mustonly contain 1 layer. Please subset the layer you would like to use for mapping.")
  }
  
  if (!stringr::str_detect(names(sraster),"strata")) {
    stop("A layer name containing 'strata' does not exist within 'sraster'.")
  }
  
  if (!stringr::str_detect(names(sraster2), "strata")) {
    stop("A layer name containing 'strata' does not exist within 'sraster2'.")
  }
  
  
  #--- map stratification rasters ---#
  
  joined <- c(sraster,sraster2)
  names(joined) <- c("strata","strata2")
  
  featuresJoin <- terra::values(joined, dataframe = TRUE)
  
  oclass <- featuresJoin %>%
    dplyr::group_by(strata, strata2) %>%
    #--- ensure NA's are transfered ---#
    dplyr::mutate(stratamapped = ifelse(is.na(strata) | is.na(strata2), NA, glue::glue('{strata}{strata2}')))
  
  #--- create lookUp table ---#
  
  lookUp <- dplyr::distinct(oclass) %>% 
    na.omit() %>%
    as.data.frame()
  
  #--- set newly stratified values ---#
  
  rout <- terra::setValues(sraster, oclass$stratamapped)
  names(rout) <- "strata"
  
  
  if(isTRUE(stack)){
    
    message("Stacking sraster, sraster2, and their combination (stratamapped).")
    
    #--- stack 3 rasters if requested ---#
    
    routstack <- c(sraster,sraster2,rout)
    names(routstack) <- c("strata","strata2","stratamapped")
    
  }
  
  #--- if not stacking rename for output ---#
  
  if(exists("routstack")){
    
    rout <- routstack
    
  }
  
  
  #--- plot if requested
  
  if(isTRUE(plot)){
    
    terra::plot(rout)
    
  }
  
  #--- write file to disc ---#
  
  if (!is.null(filename)) {
    
    #--- write file to disc depending on whether 'stack' was specified ---#
    
    if(isTRUE(stack)){
      
      terra::writeRaster(routstack, filename, overwrite = overwrite, ...)
      
    } else {
      
      terra::writeRaster(rout, filename, overwrite = overwrite, ...)
      
    }

  }
  
  #--- output details if desired ---#
  
  if (isTRUE(details)) {
    
    #--- output metrics details along with stratification raster ---#
    
    output <- list(outRaster = rout, lookUp = lookUp)
    
    #--- output samples dataframe ---#
    
    return(output)
    
  } else {
    
    #--- just output raster ---#
    
    return(rout)
  }
}