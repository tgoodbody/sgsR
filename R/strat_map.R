#' Map a raster stack of a list of rasters
#'
#' @description Map stratified rasters to a combined stratification.
#'
#' @family stratify functions
#'
#' @inheritParams strat_breaks
#' @inheritParams strat_poly
#' @param sraster spatRaster or list. Stratification raster stack or list of rasters. If \code{sraster}
#' is of class \code{list}, then it is internally converted into a raster stack.
#' @param stack Logical. Default = \code{FALSE}. If \code{TRUE}, inputs and output will be stacked:
#'  \code{strata_1, strata_2, ..., strata}.
#' @param details Logical. If \code{FALSE} (default) output is a mapped stratified spatRaster object.
#' If \code{TRUE} return a list where \code{$outRaster} is the mapped stratified raster, and
#' \code{$lookUp} is the lookup table for the stratification.
#'
#' @section Mapping:
#' The mapping algorithm will take the stratification from \code{sraster} and combine it with
#' overlying strata values across all layers. This will result in a \code{strata} attribute
#' where the values from all inputs are combined.
#'
#' i.e.
#'
#' If \code{strata_1 = 1} and \code{strata_2 = 1} then \code{strata = 11}.
#'
#' If \code{strata_1 = 2} and \code{strata_2 = 14} then \code{strata = 214}.
#' 
#' If \code{strata_1 = "A"} and \code{strata_2 = 14} then \code{strata = "A14"}.
#'
#' @examples
#' #--- load input metrics rasters ---#
#' raster <- system.file("extdata", "sraster.tif", package = "sgsR")
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
#' srasterfri <- strat_poly(
#'   poly = fri,
#'   attribute = attribute,
#'   features = features,
#'   raster = sraster
#' )
#' 
#' #--- map srasters with raster stack ---#
#' stack <- c(srasterfri, sraster)
#' strat_map(
#'   sraster = stack
#' )
#'
#' #--- map sraster with list of rasters ---#
#' rast_list <- list(srasterfri, sraster)
#' strat_map(
#'   sraster = rast_list,
#'   stack = TRUE,
#'   details = TRUE
#' )
#' @importFrom methods is
#'
#' @return A spatRaster object.
#'
#' @author Tristan R.H. Goodbody, Tommaso Trotto, Robert Hijmans
#'
#' @export

strat_map <- function(sraster,
                      stack = FALSE,
                      filename = NULL,
                      overwrite = FALSE,
                      plot = FALSE,
                      details = FALSE
) 
{
  #--- global variables ---#
  strata <- value <- NULL
  
  #--- error handling ---#
  if (!inherits(sraster, "SpatRaster") & !is.list(sraster)) {
    stop("'sraster' must be type SpatRaster or a list.", call. = FALSE)
  }

  if (is.list(sraster)) {
    if (length(sraster) <= 1) {
      stop("List must have at least 2 'SpatRaster' objects.")
    }
  }

  if (!is.logical(stack)) {
    stop("'stack' must be type logical.", call. = FALSE)
  }
  
  if (!is.logical(overwrite)) {
    stop("'overwrite' must be type logical.", call. = FALSE)
  }
  
  if (!is.logical(plot)) {
    stop("'plot' must be type logical.", call. = FALSE)
  }
  
  if (!is.logical(details)) {
    stop("'details' must be type logical.", call. = FALSE)
  }
  
  #--- list convertion into SpatRaster ---#
  if (is.list(sraster)) {
    sraster <- terra::rast(sraster)
  }
      
  #--- error handling for raster inputs ---#
  if (isFALSE(sraster@ptr$hasValues)) {
    stop("'sraster' has no values.", call. = FALSE)
  }
  
  #--- define number of input layers ---#
  nlayer <- terra::nlyr(sraster)
  
  if (nlayer <= 1) {
    stop("'sraster' must contain at least 2 layers. Please provide a 'SpatRaster' stack or a list of 'SpatRaster' objects.", call. = FALSE)
  }
  
  if (!any(grepl("strata", names(sraster)))) {
    stop("A layer name containing 'strata' does not exist within 'sraster'.", call. = FALSE)
  }
  
  #--- map stratification rasters ---#
  names(sraster) <- paste0("strata_", seq(1, nlayer))

  featuresJoin <- terra::values(sraster, dataframe = TRUE)
  
  #--- Determine index of each cell so to map values correctly without NA ---#
  
  idx <- which(complete.cases(featuresJoin))
  
  #--- check that all sraster values are the same class ---#
  classes <- c("numeric", "integer", "factor", "character")
  
  # vectorized boolean
  is_class <- sapply(seq_along(classes), function(i) {
    sapply(sapply(featuresJoin, class), identical, classes[i])
  })
  
  if (isFALSE(all(rowSums(is_class) > 0L))) {
    stop("'SpatRaster' layers must be of class 'numeric', 'integer', 'factor', or 'character'.", call. = FALSE)
  }
  
  #--- vectorize stratification across raster layers ---#
  oclass <- featuresJoin[idx,] %>%
    tidyr::unite("strata", remove = FALSE, sep = "") %>%
    dplyr::relocate(strata, .after = dplyr::last_col())

  #--- assign class to strata based on current class ---#
  if (sum(is_class[,3]) >= 1L || sum(is_class[,4]) >= 1L) { # represents classes 'factor' or 'character'
    
    oclass$strata <- as.character(oclass$strata)
    
  } else if (sum(is_class[,3]) >= 1L & sum(is_class[,4]) >= 1L) {
    
    oclass$strata <- as.character(oclass$strata)
    
  } else {
    
    oclass$strata <- as.integer(oclass$strata)
    
  }
       
  #--- create lookUp table ---#
  lookUp <- dplyr::distinct(oclass) %>%
    stats::na.omit() %>%
    as.data.frame() %>%
    dplyr::arrange(., strata)
  
  #--- make new raster with stratified values from template ---#
  odf <- matrix(nrow = nrow(featuresJoin), ncol = 1)
  
  odf[,1][idx] <- oclass$strata

  rout <- terra::setValues(sraster[[1]], odf)
  names(rout) <- "strata"
  
  #--- make raster stack with original sraster ---#
  if (isTRUE(stack)) {
    message("Stacking srasters and their combination (strata).")
    
    #--- stack rasters if requested ---#
    routstack <- c(sraster,rout)
    names(routstack) <- names(oclass)
  }
  
  #--- if stacking rename for output ---##
  if (exists("routstack")) {
    rout <- routstack
  }
  
  #--- plot if requested ---#
  if (isTRUE(plot)) {
    terra::plot(rout)
  }
  
  #--- write file to disc depending on whether 'stack' was specified ---#
  if (!is.null(filename)) {
    
    if (!is.character(filename)) {
      stop("'filename' must be type character.", call. = FALSE)
    }
    
    #--- write file to disc depending on whether 'stack' was specified ---#
    if (isTRUE(stack)) {
      terra::writeRaster(x = routstack, filename = filename, overwrite = overwrite)
      message("Output stack written to disc.")
    } else {
      terra::writeRaster(x = rout, filename = filename, overwrite = overwrite)
      message("Output raster written to disc.")
    }
  }
  
  #--- output details if desired ---#
  if (isTRUE(details)) {
    
    #--- output metrics details along with stratification raster ---#
    output <- list(raster = rout, lookUp = lookUp)
    
    #--- output samples dataframe ---#
    return(output)
  } 
  else {
    
    #--- just output raster ---#
    return(rout)
  }
}
