#' Stratify metrics raster using user defined breaks
#' @family stratify functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams strat_metrics
#' @param breaks Numeric. Vector of breakpoints for \code{metric}
#' @param breaks2 Numeric. Vector of breakpoints for \code{metric2} (if provided)
#' 
#' @return output stratification \code{spatRaster}
#' 
#' @export

strat_breaks <- function(mraster,
                      metric = NULL,
                      metric2 = NULL,
                      breaks,
                      breaks2 = NULL,
                      plot = FALSE,
                      details = FALSE)
  
{
  #--- Error management ---#
  
  if (!inherits(mraster,"SpatRaster"))
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  
  if(! metric %in% names(mraster))
    stop(paste0("mraster does not have a variable named ",metric))
  
  if (!is.numeric(breaks))
    stop("'breaks' must be type numeric")
  
  if (!is.numeric(n))
    stop("'n' must be type numeric")

  if (!is.logical(plot))
    stop("'plot' must be type logical")
  
  if (!is.logical(details))
    stop("'details' must be type logical")
  
  #--- if there is only 1 metric in the raster use it as default ---#
  
  if (terra::nlyr(mraster) == 1){
    
    rastermetric <- mraster  
    
  } else {
    
    if (is.null(metric))
      stop(" multiple layers detected in 'mraster'. Please define a 'metric' to stratify")
    
    if (!is.character(metric))
      stop("'metric' must be type character")
    
    if (any(! metric %in% names(mraster)) )
      stop(paste0("'mraster' must have an attribute named ", metric))
    
    #--- extract mraster metric ---#
    
    rastermetric <- terra::subset(mraster,metric)
    
  }
  
  #--- apply breaks to primary metric ---#
  
  #--- reclassify values based on breaks ---#
  
  breaks_m <- data.frame(from=c(-Inf,breaks,Inf)) %>%
    mutate(to = lead(from),
           becomes = seq(1:length(from))) %>%
    na.omit() %>%
    as.matrix()
  
  rcl <- terra::classify(rastermetric,breaks_m)
  names(rcl) <- "strata"
  
  
  #--- if secondary metric is provided ---#
  
  if (!is.null(metric2)){
    
    if (!is.character(metric2))
      stop("'metric2' must be type character")
    
    if (is.null(breaks2))
      stop("If using metrics2 to stratify, 'breaks2' must be defined")
    
    if (any(! metric2 %in% names(mraster)) )
      stop(paste0("'mraster' must have an attribute named ", metric2))
    
    #--- extract metric2 from mraster ---#
    
    rastermetric2 <- terra::subset(mraster, metric2)
    
    #--- reclassify values based on breaks ---#
    
    breaks2_m <- data.frame(from=c(-Inf,breaks2,Inf)) %>%
      mutate(to = lead(from),
             becomes = seq(1:length(from))) %>%
      na.omit() %>%
      as.matrix()
    
    rcl2 <- terra::classify(rastermetric2,breaks2_m)
    names(rcl2) <- "strata2"
    
    #--- stack rcl and rcl2 
    
    rstack <- c(rcl, rcl2)
    
    breaks_c <- terra::as.data.frame(rstack, xy=TRUE)
    
    breaks_c <- breaks_c %>%
      group_by(strata,strata2) %>%
      #--- establish newly formed unique class ---#
      mutate(class = cur_group_id()) %>%
      ungroup()
    
    idx <- cellFromXY(rcl,cbind(breaks_c$x,breaks_c$y))
    
    #--- convert back to original mraster extent ---#
    rcl[idx] <- breaks_c$class
    
  }
  
  
  if (isTRUE(plot)){
    
    data <- terra::as.data.frame(rastermetric)
    names(data) <- "metric"
    
    #--- plot histogram of metric with associated break lines ---#
    
    p1 <- ggplot(data,aes(metric)) +
      geom_histogram() +
      geom_vline(xintercept = breaks, linetype = "dashed") +
      ggtitle(paste0(metric, " histogram with defined breaks"))
    
    if (!is.null(metric2)){
      
      data2 <- terra::as.data.frame(rastermetric2)
      names(data2) <- "metric2"
    
    p2 <- ggplot(data2,aes(metric2)) +
      geom_histogram() +
      geom_vline(xintercept = breaks2, linetype = "dashed") +
      ggtitle(paste0(metric2, " histogram with defined breaks"))
    
    suppressMessages(print(ggarrange(p1, p2, ncol = 1, nrow =2)))
    
    } else {
      
      suppressMessages(print(p1))
      
    }
    
    #--- set colour palette ---#
    
    terra::plot(rcl, main = 'User break defined strata',type="classes")
    
  }
  
  #--- Output based on 'details' to return raster alone or list with details ---#
  
  if ( isTRUE(details) ){
    
    #--- output break points raster with associated breaks ---#
    
    if (!is.null(metric2)){
    
    breaks_rcl <- list(details = list(breaks = breaks, 
                                      breaks2 = if(!missing(breaks2)) breaks2), 
                       raster = rcl)
    
    return(breaks_rcl)
    
    } else {
      
      #--- output break points raster with associated breaks ---#

        breaks_rcl <- list(details = breaks, 
                           raster = rcl)
        
        return(breaks_rcl)
        
      
    }
    
    
  } else {
    
    #--- just output raster ---#
    
    return(rcl)
    
  }
  
}


