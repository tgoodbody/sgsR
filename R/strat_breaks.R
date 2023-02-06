#' Breaks stratification
#'
#' @description Stratify metrics raster using user defined breaks
#' @family stratify functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams sample_systematic
#' @param mraster Spatraster. Raster to stratify. Layers in \code{mraster} must match the number of
#' \code{breaks} vectors provided.
#' @param breaks Numeric. Vector of breakpoints for each layer of \code{mraster}. If \code{mraster} has multiple layers,
#' \code{breaks} must be a list with an equal number of objects.
#' @param map Logical. Map individual stratified layers to a combined stratification. Will output a multi-layer
#' \code{SpatRaster} with individual stratifications for each \code{mraster} layer and an additional mapped stratification
#'  named \code{"strata"}.
#' @param filename Character. Path to write stratified raster to disc.
#' @param overwrite Logical. Specify whether \code{filename} should be overwritten on disc.
#'
#' @return Returns an output stratification \code{spatRaster} or a list when \code{details = TRUE}.
#'
#' When a list is returned:
#' \enumerate{
#' \item \code{raster} is a stratified \code{spatRaster} based on quantiles. If \code{stack = TRUE} will
#' be the number of layers of \code{mraster} plus the final output
#' \item \code{breaks} is a list output of \code{breaks}
#' \item \code{plot} is a \code{ggplot} histogram object showing distribution(s) and break point(s).
#' }
#'
#' @examples
#' #--- Load raster ---#
#' r <- system.file("extdata", "mraster.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' #--- create vector breaks ---#
#' br.zq90 <- c(3, 5, 11, 18)
#' br.pz2 <- c(20, 40, 60, 80)
#'
#' strat_breaks(
#'   mraster = mr$zq90,
#'   breaks = br.zq90
#' )
#' 
#'
#' strat_breaks(
#'   mraster = mr[[1:2]],
#'   breaks = list(br.zq90, br.pz2),
#'   details = TRUE
#' )
#' @author Tristan R.H. Goodbody
#'
#' @export

strat_breaks <- function(mraster,
                         breaks,
                         map = FALSE,
                         plot = FALSE,
                         details = FALSE,
                         filename = NULL,
                         overwrite = FALSE
) {
    
  #--- Set global vars ---#
  strata <- x <- y <- value <- val <- NULL
  
  #--- Error management ---#
    
  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster.", call. = FALSE)
  }
    
  if (!all(sapply(breaks, is.numeric))) {
    stop("'breaks' must be type numeric.", call. = FALSE)
  }
  
  if (!is.logical(map)) {
    stop("'map' must be type logical.", call. = FALSE)
  }
  
  if (!is.logical(plot)) {
    stop("'plot' must be type logical.", call. = FALSE)
  }
  
  if (!is.logical(details)) {
    stop("'details' must be type logical.", call. = FALSE)
  }
  
  #--- get number of layers ---#
  nlayer <- terra::nlyr(mraster)
  
  #--- if there is only 1 band in mraster use it as default ---#
  if(nlayer > 1){
    
    if(!inherits(breaks, "list")){
      stop("`breaks` must be a list of numeric vectors of the same length as `mraster`.", call. = FALSE)
    }
    
    if (nlayer != length(breaks)) {
      stop("`mraster` and `breaks` must have the same number of layers & objects.", call. = FALSE)
    }
    
  } else{
    
    if(!is.list(breaks)){
      breaks <- list(breaks)
    }
    
  }
  
  #--- check that breaks values are not < or > raster min/max raster values ---#
  
  minmax <- terra::minmax(mraster)
  
  mins <- which(unlist(mapply(function(x,y) x > y, breaks, minmax[1,])) == FALSE)
  maxs <- which(unlist(mapply(function(x,y) x < y, breaks, minmax[2,])) == FALSE)
  
  if (length(mins) != 0) {
    stop("'breaks' contains values < the minimum corresponding 'mraster' value.", call. = FALSE)
  }
  
  if (length(maxs) != 0) {
    stop("'breaks' contains values > the maximum corresponding 'mraster' value.", call. = FALSE)
  }
  
  #--- convert mraster to list ---#
  mrl <- as.list(mraster)
  
  #--- vectorize vect_breaks to determine raster break points ---#
  rstack <- terra::rast(mapply(calculate_breaks, mraster = mrl, breaks = breaks))
  
  #--- rename to append raster metric name ---#
  names(rstack) <- paste0("strata_",names(mraster))
  
  if(nlayer > 1){
    #--- establish newly formed unique strata ---#
    breaks_c <- terra::as.data.frame(rstack, xy = TRUE)%>%
      dplyr::group_by(dplyr::across(tidyr::starts_with("strata"))) %>%
      dplyr::mutate(strata = dplyr::cur_group_id()) %>%
      dplyr::ungroup()
    
    lookUp <- dplyr::arrange(breaks_c, strata) %>% 
      dplyr::select(-x,-y) %>% 
      unique()
    
    idx <- terra::cellFromXY(mraster[[1]], cbind(breaks_c$x, breaks_c$y))
    
    rcl <- mraster[[1]]
    #--- convert back to original extent ---#
    rcl[idx] <- breaks_c$strata
  
    names(rcl) <- "strata"
    
    if(isTRUE(map)){
      
      message("Mapping stratifications.")
      
      rcl <- c(rstack,rcl)
      
    } else {
      
      rcl <- rstack
      
      names(rcl) <- rep("strata",nlayer)
      
    }
    
  } else {
    
    rcl <- rstack
    
    names(rcl) <- "strata"
  }
  
  if (isTRUE(plot)) {
    
    brs <- mapply(function(x,y) data.frame(var = x, brk = y), names(mraster), breaks)
    
    brs <- do.call("rbind",apply(X = brs, MARGIN = 2, FUN = function(x) data.frame(names = x$var, val = x$brk)))
    
    p <- terra::as.data.frame(mraster) %>%
      tidyr::pivot_longer(dplyr::everything(),names_to = "names") %>%
      ggplot2::ggplot(ggplot2::aes(value)) +
      ggplot2::geom_histogram() +
      ggplot2::geom_vline(linetype = "dashed", data = brs, mapping = ggplot2::aes(xintercept = val)) +
      ggplot2::facet_wrap(~names, scales = "free")
    
    suppressMessages(print(p))
    
    #--- set colour palette ---#
    
    terra::plot(rcl)
  }
  
  #--- write file to disc ---#
  
  if (!is.null(filename)) {
    terra::writeRaster(x = rcl, filename = filename, overwrite = overwrite)
    message("Output raster written to disc.")
  }
  
  #--- Output based on 'details' to return raster alone or list with details ---#
  
  if (isTRUE(details)) {

    #--- output break points raster with associated breaks ---#
    breaks_rcl <- list(
      breaks = if (exists("lookUp")) lookUp,
      raster = rcl,
      plot = if (exists("p")) p
    )
    
    return(breaks_rcl)
  } else {
    
    #--- just output raster ---#
    
    return(rcl)
  }
}

