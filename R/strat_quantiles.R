#' Quantiles stratification
#'
#' @description Stratify metric raster using metric quantiles.
#' @family stratify functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams strat_breaks
#' @param nStrata Numeric. A positive integer representing the number of desired classes 
#' or a numeric vector of probabilities with values between \code{0-1}. If \code{mraster} has multiple layers,
#' \code{nStrata} must be a list with an equal number of objects. 
#'
#' @return Returns an output stratification \code{spatRaster} or a list when \code{details = TRUE}.
#'
#' When a list is returned:
#' \enumerate{
#' \item \code{details} lookUp table for stratification(s).
#' \item \code{raster} is a stratified \code{spatRaster} based on quantiles
#' \item \code{plot} is a \code{ggplot} histogram / scatter plot object (depends on whether metric2 was supplied).
#' Histogram shows distribution and break points.
#' }
#'
#' @examples
#' #--- Load raster and existing plots---#
#' r <- system.file("extdata", "mraster.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' strat_quantiles(
#'   mraster = mr$zq90,
#'   nStrata = 4
#' )
#'
#' strat_quantiles(
#'   mraster = mr[[1:2]],
#'   nStrata = list(c(0.2,0.4,0.8), 3),
#'   map = TRUE
#' )
#'   
#' @author Tristan R.H. Goodbody
#'
#' @export

strat_quantiles <- function(mraster,
                            nStrata,
                            map = FALSE,
                            plot = FALSE,
                            details = FALSE,
                            filename = NULL,
                            overwrite = FALSE) {

  #--- Set global vars ---#
  strata <- x <- y <- value <- val <- NULL

  #--- error handling ---#
  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster.", call. = FALSE)
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
    
    if (nlayer != length(nStrata)) {
      stop("`mraster` and `nStrata` must have the same number of layers & objects.", call. = FALSE)
    }
    
    if(!inherits(nStrata, "list")){
      stop("`nStrata` must be a list of numeric vectors of the same length as `mraster`.", call. = FALSE)
    }
    
  } else{
    
    if(!is.list(nStrata)){
      nStrata <- list(nStrata)
    }
    
  }
  
  #--- convert mraster to list ---#
  mrl <- as.list(mraster)
  
  #--- vectorize vect_breaks to determine raster break points ---#
  rastbreaks <- mapply(calculate_quantile_breaks, mraster = mrl, nStrata = nStrata)
  
  #--- take all raster outputs from vectorization ---#
  rstack <- terra::rast(rastbreaks[c(TRUE, FALSE)])
  
  #--- rename to append raster metric name ---#
  names(rstack) <- paste0("strata_",names(mraster))
  
  if(nlayer > 1){
    #--- establish newly formed unique class ---#
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
    
    #--- take numeric breaks for plotting ---#
    brs <- do.call("rbind",rastbreaks[c(FALSE, TRUE)])

    brs <- brs[is.finite(brs$val),]

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
    #--- output metrics details along with stratification raster ---#

    if (isTRUE(plot)){
    out <- list(
      details = lookUp,
      raster = rcl,
      plot = p
    )
    } else {
      out <- list(
        details = lookUp,
        raster = rcl
      )
    }

    return(out)
  } else {

    #--- just output raster ---#

    return(rcl)
  }
}
