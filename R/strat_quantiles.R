#' Quantiles stratification
#' 
#' @description Stratify metric raster using metric quantiles.
#' @family stratify functions
#'
#' @inheritParams strat_kmeans
#' @param metric Character. Name of primary metric to stratify. If 
#' \code{mraster} is has 1 layer it is taken as default.
#' @param metric2 Character. Name of secondary metric to stratify.
#' @param nstrata2 Numeric.  Number of secondary strata within \code{nstrata}.
#' @param samp Numeric. For plotting - Determines proportion of cells 
#' for strata visualization. Lower values reduce processing time.
#'
#' @importFrom magrittr %>%
#' @importFrom methods is
#' 
#' @return output stratification \code{spatRaster}, or a list when \code{details = TRUE}.
#' 
#' @export

strat_quantiles <- function(mraster,
                          metric = NULL,
                          metric2 = NULL,
                          nstrata,
                          nstrata2 = NULL,
                          plot = FALSE,
                          samp = 1,
                          details = FALSE){
  
  #--- Set global vars ---#
  
  class1 <- class2 <- NULL
  
  #--- error handling ---#
  if (!inherits(mraster,"SpatRaster"))
    stop("all specified bands must be type SpatRaster", call. = FALSE)
  
  if (!is.numeric(nstrata))
    stop("'nstrata' must be type numeric")
  
  if (!is.logical(plot))
    stop("'plot' must be type logical")
  
  if (!is.numeric(samp))
    stop("'samp' must be type numeric")
  
  if (!is.logical(details))
    stop("'details' must be type logical")
  
  if (is.null(metric2)){
    if (!is.null(nstrata2))
      message("You are stratifying with only 1 metric but specified 'nstrata2' - ignoring.")
    
    #--- if there is only 1 metric in the raster use it as default ---#
    
    if (terra::nlyr(mraster) == 1){
      
      #--- Extract values from mraster ---#
      
      vals <- terra::values(mraster)  
      
      #--- set name of raster band to 'metric' ---#
      
      metric <- names(mraster)
      
    } else {
      
      if (is.null(metric))
        stop(" multiple layers detected in 'mraster'. Please define a 'metric' to stratify")
      
      if (!is.character(metric))
        stop("'metric' must be type character")
      
      if (any(! metric %in% names(mraster)) )
        stop(paste0("'mraster' must have an attribute named ", metric))
      
      #--- Extract values from mraster ---#
      
      vals <- terra::subset(mraster,metric) %>%
        terra::values()
      
    }
    
    vals[!is.finite(vals)] <- NA
    
    #--- Determine index of each cell so to map values correctly without NA ---#
    
    idx <- !is.na(vals)
    
    #--- Remove NA / NaN / Inf values ---#
    
    df <- vals %>%
      as.data.frame() %>%
      dplyr::filter(!is.na(.))
    
    metric <- ggplot2::ensym(metric)
    
    #--- Split metric distribution in to number specified by 'breaks' ---#
    
    dfc <- df %>%
      dplyr::mutate(class = dplyr::ntile(!!metric,nstrata))
    
    #--- convert back to original mraster extent ---#
    
    vals[idx] <- dfc$class
    
    #--- set newly stratified values ---#
    
    rout <- terra::setValues(mraster[[1]],vals)
    names(rout) <- "strata"
    
    
  }
  
  if (!is.null(metric2)){
    
    if (!is.character(metric2))
      stop("'metric2' must be type character")
    
    if (is.null(nstrata2))
      stop("If using 2 metrics to stratify, 'nstrata2' must be defined")
    
    if (any(! metric2 %in% names(mraster)) )
      stop(paste0("'mraster' must have an attribute named ", metric2))
    
    #--- Extract values from mraster ---#
    
    vals <- terra::subset(mraster,c(metric,metric2)) %>%
      terra::values()
    
    vals[!is.finite(vals)] <- NA
    
    #--- Determine index of each cell so to map values correctly without NA ---#
    
    idx <- is.finite(vals[,1]) & is.finite(vals[,2])
    
    #--- Remove NA / NaN / Inf values ---#
    
    df <- vals %>%
      as.data.frame() %>%
      dplyr::filter(!is.na(.))
    
    metric <- ggplot2::ensym(metric)
    metric2 <- ggplot2::ensym(metric2)
    
    #--- Split metric distribution in to number specified by 'breaks' ---#
    
    dfc <- df %>%
      #--- define nstrata classes ---#
      dplyr::mutate(class1 = dplyr::ntile(!!metric,nstrata)) %>%
      #--- group by class to sub stratify ---#
      dplyr::group_by(class1) %>%
      #--- define nstrata2 classes ---#
      dplyr::mutate(class2 = dplyr::ntile(!!metric2,nstrata2)) %>%
      #--- combine classes ---#
      dplyr::group_by(class1,class2) %>%
      #--- establish newly formed unique class ---#
      dplyr::mutate(class = dplyr::cur_group_id())
    
    #--- convert back to original mraster extent ---#
    
    vals[,1][idx] <- dfc$class
    
    #--- set newly stratified values ---#
    
    rout <- terra::setValues(mraster[[1]],vals[,1])
    names(rout) <- "strata"
    
  }
    
  if (isTRUE(plot)){
    if (samp > 1 | samp < 0)
      stop("'samp' must be between 0 and 1")
    
    if(is.null(metric2)){
      
      #--- output histogram with quantile breaks ---#
      
      #--- determine numeric breaks ---#
      
      breaks <- dfc %>% 
        dplyr::group_by(class) %>%
        dplyr::summarize(breaks = max(!!metric)) %>%
        dplyr::select(breaks) %>%
        as.data.frame()
      
      breaks <- breaks[1:(nstrata-1),]
      
      df.p <- dfc %>% 
        dplyr::select(metric)
      
      #--- plot histogram of metric with associated break lines ---#
      
      p1 <- ggplot2::ggplot(df.p,ggplot2::aes(!!metric)) +
        ggplot2::geom_histogram() +
        ggplot2::geom_vline(xintercept = breaks, linetype = "dashed") +
        ggplot2::ggtitle(paste0(metric, " histogram with defined breaks"))
      
      print(p1)
      
    } else {
    
    #--- set up colour palette ---#
    
    ncol <- nstrata * nstrata2
    qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'seq',]
    col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    terra::plot(rout, main = 'Classes', col=sample(col_vector, ncol))
    
    coordsgrps <- dfc %>%
      dplyr::group_by(class) %>%
      dplyr::arrange(class) %>%
      tidyr::nest() %>%
      dplyr::ungroup()
    
    q <- classPlot(dfc = dfc,
                   coordsgrps = coordsgrps,
                   metric = metric,
                   metric2 = metric2,
                   samp = samp)
    
    print(q)
    
  }
    
    #--- Output based on 'details' to return raster alone or list with details ---#
    
    if ( isTRUE(details) ){
      
      #--- output metrics details along with stratification raster ---#
      
      out <- list(details = dfc, 
                  raster = rout)
      
      return(out)
      
      
    } else {
      
      #--- just output raster ---#
      
      return(rout)
      
    }
    
  }
  
  return(rout)
  
}


