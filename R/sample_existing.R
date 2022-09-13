#' Sample existing
#' 
#' @description Sub-sample an existing sample using \code{\link[clhs]{clhs}} functionality.
#'
#' @family sample functions
#' 
#' @inheritParams sample_systematic
#' @inheritParams extract_strata
#' @inheritParams sample_clhs
#' 
#' @param raster spatRaster. Raster used to define population distributions.
#' @param ... Additional arguments for clhs sampling. See \code{\link[clhs]{clhs}}.
#' 
#' @return An sf object of samples
#' 
#' @note 
#' 
#' If providing only \code{existing} - all attributes will be used for sampling. Remove attributes not indented for sampling
#' prior to using this algorithm.
#' 
#' @examples 
#' #--- Load raster ---#
#' r <- system.file("extdata", "mraster.tif", package = "sgsR")
#' mr <- terra::rast(r)
#' 
#' #--- generate an existing sample ---#
#' e <- sample_systematic(raster = mr, cellsize = 200)
#' 
#' #--- perform sub-sampling ---#
#' sample_existing(existing = e,
#'                 raster = mr,
#'                 nSamp = 50)
#'                 
#' #--- extract metrics to sample ---#
#' e <- extract_metrics(mr, e)
#' 
#' #--- perform sub-sampling ---#
#' sample_existing(existing = e,
#'                 nSamp = 30, 
#'                 plot = TRUE)
#' 
#' @author Tristan R.H. Goodbody
#' 
#' @export

sample_existing <- function(existing,
                            nSamp,
                            raster = NULL,
                            cost = NULL,
                            access = NULL,
                            buff_inner = NULL,
                            buff_outer = NULL,
                            iter = 10000,
                            plot = FALSE,
                            details = FALSE,
                            filename = NULL,
                            overwrite = FALSE,
                            ...
){
  
  #--- set global variables ---#
  
  x <- y <- X <- Y <- type <- value <- NULL
  
  #--- error handling ---#
  
  if (!inherits(existing, "data.frame") && !inherits(existing, "sf")) {
    stop("'existing' must be a data.frame or sf object.", call. = FALSE)
  }
  
  if (!is.null(raster) & !inherits(raster, "SpatRaster")) {
    stop("'raster' must be type SpatRaster.", call. = FALSE)
  }
  
  if (!is.numeric(nSamp)) {
    stop("'nSamp' must be type numeric of type integer or double.", call. = FALSE)
  }
  
  if (nSamp >= nrow(existing)){
    stop("'nSamp' must be less than the total number of 'existing' samples.", call. = FALSE)
  }
  
  if (!is.numeric(iter)) {
    stop("'iter' must be type numeric.", call. = FALSE)
  }
  
  if (iter <= 0) {
    stop("'iter' must be  >= 0.", call. = FALSE)
  }
  
  if (!is.logical(plot)) {
    stop("'plot' must be type logical.", call. = FALSE)
  }
  
  if (!is.logical(details)) {
    stop("'details' must be type logical.", call. = FALSE)
  }
  
  #--- Prepare existing sample data ---#
  if(!inherits(existing, "sf")){
    
    if (any(!c("X", "Y") %in% colnames(existing))) {
      
      #--- if coordinate column names are lowercase change them to uppercase to match requirements ---#
      
      if (any(c("x", "y") %in% colnames(existing))) {
        existing <- existing %>%
          dplyr::rename(
            X = x,
            Y = y
          )
        
        message("'existing' column coordinate names are lowercase - converting to uppercase.")
      } else {
        
        #--- if no x/y columns are present stop ---#
        
        stop("'existing' must have columns named 'X' and 'Y'.", call. = FALSE)
      }
    }
    
    #--- if raster is supplied
    if(!is.null(raster)){
      
      #--- check to see if 'existing' does not contain attributes with the same names as 'raster' ---#
      if(!all(names(raster) %in% names(existing))){
        
        message("'existing' does not contain attributes with the same names as 'raster'. Extracting metrics.")
        
        existing <- existing %>%
          sf::st_as_sf(., coords = c("X", "Y"))
        
        #--- assign sraster crs to spatial points object ---#
        sf::st_crs(existing) <- crs
        
        #--- access constraint ---#
        
        if(!is.null(access)){
          
          masked <- mask_existing(access = access, existing = existing, buff_inner = buff_inner, buff_outer = buff_outer)
          
          existing <- masked$samples
          
        }
        
        #--- extract covariates at existing sample locations ---#
        existingdf <- extract_metrics(mraster = raster, existing = existing, data.frame = TRUE) %>%
          na.omit()
        
      }
      
    }

  } else {
    
    #--- check to see if 'existing' does not contain attributes with the same names as 'raster' ---#
    if(!all(names(raster) %in% names(existing))){
      
      message("'existing' does not contain attributes with the same names as 'raster'. Extracting metrics.")
      
      #--- access constraint ---#
      
      if(!is.null(access)){
        
        masked <- mask_existing(access = access, existing = existing, buff_inner = buff_inner, buff_outer = buff_outer)
        
        existing <- masked$samples
        
      }
      
      #--- extract covariates at existing sample locations ---#
      existingdf <- extract_metrics(mraster = raster, existing = existing, data.frame = TRUE) %>%
        na.omit()
      
    } else {
      
      xy <- existing %>%
        sf::st_coordinates(.)
      
      existingdf <- existing %>%
        sf::st_drop_geometry(.) %>%
        cbind(xy, .)
    }
    
  }
  
  #--- incorporate cost constraint ---#
  
  if (!is.null(cost)) {
    
    if(!is.character(cost) & !is.numeric(cost)){
      stop("'cost' must be either type numeric or character.", call. = FALSE)
    }
    
    if(is.null(raster)){ 
      
      costLoc <- existing 
      
      v <- "existing"
    
    } else {
      
      costLoc <- raster
      
      v <- "raster"
    }
    
    if (is.numeric(cost)) {
      if (cost > length(names(costLoc)) | cost < 0) {
        stop(paste0("'cost' index doest not exist within '", v, "'."), call. = FALSE)
      }
      
    } else {
      
      if(length(which(names(costLoc) == cost)) == 0){
        
        stop(paste0("No layer named '",cost,"' exists in '", v, "'."), call. = FALSE)
        
      } else {
        
        cost <- which(names(costLoc) == cost)
        
      }
    }
    
    message(paste0("Using `", names(costLoc)[cost], "` as sampling constraint."))
    
  } else {
    
    costLoc <- NULL
    
  }
  
  #--- determine if sampling will be performed using metric attributes in 'existing' only or with the addition of 'raster' metrics ---#
  
  if(is.null(raster)){
    
    #--- test if existing has attributes other than geometry ---#
    if(ncol(existing) <= 1){
      if(names(existing) == "geometry"){
        stop("'existing' has no attributes that can be sampled.", call. = FALSE)
      }
    }
    
    message("Sub-sampling based on 'existing' metric distributions.")

    all <- existing %>%
      sf::st_drop_geometry(.)
    
    #--- sampling ---#
    
    sidx <- (nrow(all)-nrow(existing)+1):nrow(all)
    
    if(isTRUE(details)){
      
      #--- output clhs information to be supplied in 'details' list output ---#
      
      clhsOut <- clhs::clhs(x = all, size = nSamp, iter = iter, cost = cost, can.include = sidx, simple = FALSE, ...)
      
      outIdx <- clhsOut$index_samples
      #--- extract sampled rows from existing ---#
      
      samples <- existing[outIdx,]

    } else {
      
      outIdx <- clhs::clhs(x = all, size = nSamp, iter = iter, cost = cost, can.include = sidx, ...)
      
      #--- extract sampled rows from existing ---#
      
      samples <- existing[outIdx,]
    }
  
    } else {
    
      message("Sub-sampling based on 'raster' distributions.")
    
      #--- determine crs of input raster ---#
      crs <- terra::crs(raster, proj = TRUE)
      
      # #--- extract covariates data from raster ---#
      vals <- terra::as.data.frame(raster, xy = TRUE, row.names = FALSE) %>%
        dplyr::rename(
          X = x,
          Y = y
        ) %>%
        stats::na.omit()
    
      #--- select the variables existing in raster for sampling ---#
      
      e <- existingdf %>%
        dplyr::select(X,Y, names(raster))
      
      all <- rbind(vals, e)
        
      
      #--- sampling ---#
      
      sidx <- (nrow(all)-nrow(e)+1):nrow(all)
      
      if(isTRUE(details)){
        
        #--- output clhs information to be supplied in 'details' list output ---#
        
        clhsOut <- clhs::clhs(x = all[,3:ncol(all)], size = nSamp, iter = iter, cost = cost, can.include = sidx, simple = FALSE, ...)
        
        outIdx <- clhsOut$index_samples
        #--- extract sampled rows from existing ---#
        
        samples <- all[outIdx,]
        
      } else {
        
        outIdx <- clhs::clhs(x = all[,3:ncol(all)], size = nSamp, iter = iter, cost = cost, can.include = sidx, ...)
        
        #--- extract sampled rows from existing ---#
        
        samples <- all[outIdx,]
      }
      
      #--- extract sampled rows from existing ---#
      
      samples <- samples %>%
        sf::st_as_sf(., coords = c("X", "Y")) %>%
        sf::st_set_crs(., sf::st_crs(existing))
    
  }
  
  #--- plot ---#
  
  if(isTRUE(plot)){
    
    #--- generate data.frames for original population and selected samples ---#
    pop <- all %>%
      dplyr::mutate(type = "population")
    
    all <- all[outIdx,] %>%
      dplyr::mutate(type = "sample") %>%
      rbind(pop, .)
    
    #--- plotting ---#
    if(is.null(raster)){
      
      pecdf <- all %>% 
        dplyr::select(-dplyr::all_of(names(costLoc)[cost])) %>%
        tidyr::pivot_longer(c(!type), names_to = "metric") %>%
        ggplot2::ggplot(ggplot2::aes(value, class = type, colour = type)) + 
        ggplot2::stat_ecdf(geom = "step") +
        ggplot2::facet_grid(.~ metric, scales = "free") +
        ggplot2::ggtitle(label = "Empirical cumulative distributions for existing samples") +
        ggplot2::xlab("Sampling metrics") +
        ggplot2::ylab("Cumulative distribution")
      
    } else {

      if(!is.null(access)){
        #--- plot sample locations and raster output ---#
        terra::plot(raster[[1]])
        suppressWarnings(plot(masked$buffer, add = TRUE, alpha = 0.1))
        suppressWarnings(plot(samples, add = TRUE, col = "black"))
      } else {
        #--- plot sample locations and raster output ---#
        terra::plot(raster[[1]])
        suppressWarnings(plot(samples, add = TRUE, col = "black"))
      }
      
      #--- generate ecdf curves comparing population and sample ---#
      pecdf <- all %>% 
        dplyr::select(-X,-Y,-dplyr::all_of(names(costLoc)[cost])) %>%
        tidyr::pivot_longer(c(!type), names_to = "metric") %>%
        ggplot2::ggplot(ggplot2::aes(value, class = type, colour = type)) + 
        ggplot2::stat_ecdf(geom = "step") +
        ggplot2::facet_grid(.~ metric, scales = "free") +
        ggplot2::ggtitle(label = "Empirical cumulative distributions for existing samples") +
        ggplot2::xlab("Sampling metrics") +
        ggplot2::ylab("Cumulative distribution")
    }
  
    print(pecdf)
  
  }
  
  #--- write to disc ---#
  
  if (!is.null(filename)) {
    
    if (!is.character(filename)) {
      stop("'filename' must be a file path character string.", call. = FALSE)
    }
    
    if (!is.logical(overwrite)) {
      stop("'overwrite' must be type logical.", call. = FALSE)
    }
    
    if (file.exists(filename) & isFALSE(overwrite)) {
      stop(paste0("'",filename, "' already exists and overwrite = FALSE."), call. = FALSE)
    }
    
    sf::st_write(samples, filename, delete_layer = overwrite)
    message("Output samples written to disc.")
  }
  
  #--- details ---#

  if(isTRUE(details)){
   
    out <- list(samples = samples, 
                population = all, 
                clhsOut = clhsOut, 
                plot = if (exists("pecdf")) pecdf)
    
    return(out)
     
  } else {
    
    return(samples)
  }
  
}
