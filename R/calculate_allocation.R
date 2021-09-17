#' Determine required samples in strata
#'
#' @details Determine how many samples to allocate within strata
#'
#' @family calculate functions
#'
#' @inheritParams sample_srs
#' @inheritParams sample_strat
#' @param allocation Character. Allocation algorithm to be used. Either \code{prop} (default) for proportional allocation
#' or \code{optim} for optimal allocation or \code{equal} for equal number of samples for each strata.
#' @param mraster spatRaster. ALS metrics raster. Required when \code{allocation = optim}.
#' @param metric Numeric/Character. Index or name of primary covariate within mraster to stratify.
#' Required when \code{allocation = optim}.
#' @param force Logical. \code{Default = FALSE} - force \code{nSamp} to be exactly the user defined value
#' in cases where nSamp and \code{sraster} strata count are not equally divisible. Has no effect when \code{existing}
#' is provided.
#'
#' @return data.frame of:
#' \itemize{
#' \item{strata} - Strata ID.
#' \item{total} - Total samples to be allocated (positive) or removed (negative).
#' \item{need} - Total required samples per strata.
#' }
#' 
#' @section Details:
#' \itemize{
#' \item{nSamp = Desired sample count}
#' \item{Nh = Count of cells within strata 'h'}
#' \item{N = Total cound of cells}
#' \item{\eqn{\sigma}h = Sample standard deviation}
#' \item{nh = samples allocated with strata 'h'}
#' }
#' 
#' \eqn{nh = nSamp * Nh / N}
#' 
#' \eqn{nh = nSamp  Nh * Nh * \sigma h / \sum Nk * \sigma k}
#' 
#'
#' @examples
#' #--- Load strata raster and existing samples---#
#' r <- system.file("extdata", "kmeans.tif", package = "sgsR")
#' sr <- terra::rast(r)
#'
#' e <- system.file("extdata", "existing.shp", package = "sgsR")
#' e <- sf::st_read(e)
#'
#' #--- perform grid sampling ---#
#' calculate_allocation(sraster = sr, 
#'                      nSamp = 200)
#'
#' calculate_allocation(sraster = sr, 
#'                      nSamp = 200,
#'                      force = TRUE)
#'                      
#' #--- extract strata from existing samples ---#
#' e.sr <- extract_strata(sraster = sr, 
#'                        existing = e)
#' 
#' calculate_allocation(sraster = sr, 
#'                      nSamp = 200, 
#'                      existing = e.sr)
#' 
#' #--- Load mraster for optimal allocation ---#                     
#' mr <- system.file("extdata", "wall_metrics_small.tif", package = "sgsR")
#' mr <- terra::rast(mr)
#'                      
#' calculate_allocation(sraster = sr, 
#'                      nSamp = 200, 
#'                      existing = e.sr,
#'                      allocation = "optim",
#'                      mraster = mr,
#'                      metric = 1,
#'                      force = TRUE)
#'                      
#'                      
#' @author Tristan R.H. Goodbody
#'
#' @export


calculate_allocation <- function(sraster,
                                 nSamp,
                                 allocation = "prop",
                                 mraster = NULL,
                                 metric = NULL,
                                 existing = NULL,
                                 force = FALSE) {
  
  #--- Error management ---#
  if (!inherits(sraster, "SpatRaster")) {
    stop("'sraster' must be type SpatRaster", call. = FALSE)
  }
  
  if (!is.numeric(nSamp)) {
    stop("'nSamp' must be type numeric")
  }
  
  if (allocation != "prop" && allocation != "optim" && allocation != "equal") {
    stop("Unknown allocation: '", allocation, "' selected. Please use 'prop' (default) or 'optim' or 'equal'")
  }
  
  if (!is.logical(force)) {
    stop("'force' must be type logical")
  }
  
  #--- set global vars ---#
  
  strata <- count <- freq <- total <- eTotal <- NULL
  
  #--- determine crs of input sraster ---#
  crs <- terra::crs(sraster, proj=TRUE)
  
  #--- determine which allocation algorithm to use ---#
  
  #--- proportional allocation ---#
  
  if(allocation == "prop"){
    
    message("Implementing porportional allocation of samples")
    
    if(!is.null(mraster)){
      message("'allocation = prop' - ignoring 'mraster'")
    }
    
    if(!is.null(metric)){
      message("'allocation = prop' - ignoring 'metric'")
    }
    
    vals <- terra::values(sraster) %>%
      as.data.frame()
    
    names(vals) <- "strata"
    
    #--- determine number of samples within each strata ---#
    toSample <- vals %>%
      stats::na.omit() %>%
      dplyr::group_by(strata) %>%
      dplyr::summarize(count = dplyr::n()) %>%
      dplyr::mutate(
        freq = count / sum(count),
        total = freq * nSamp
      ) %>%
      #--- if a value equates to <1 it will have 0 samples --- change 0 to 1 ---#
      
      #########################################
    #### What other method could be used ####
    #########################################
    
    dplyr::mutate(total = replace(total, total < 1, 1)) %>%
      dplyr::mutate(total = round(total)) %>%
      dplyr::select(strata, total) %>%
      as.data.frame()
    
  }
  
  #--- optimal allocation ---#
  
  if(allocation == "optim"){
    
    #--- error handling when allocation algorithm is 'optim' ---#
    
    if(is.null(mraster)){
      stop("'mraster' must be supplied if 'allocation = optim'.")
    }
    
    if (!inherits(mraster, "SpatRaster")) {
      stop("'mraster' must be type SpatRaster", call. = FALSE)
    }
    
    #--- determine numeric or character index of metric for variability calculation ---#
    
    if (terra::nlyr(mraster) == 1) {
      
      #--- set name of raster band to 'metric' ---#
      
      metric <- names(mraster)
    } else {
      
      #--- subset metric based on whether it is a character of number ---#
      
      if (is.null(metric)) {
        stop(" multiple layers detected in 'mraster'. Define 'metric' to calculate stratum variability.")
      } else {
        
        #--- Numeric ---#
        
        if (is.numeric(metric)) {
          if ((metric) > (terra::nlyr(mraster)) | metric < 0) {
            stop("'metric' index doest not exist within 'mraster'")
          }
          
          metric <- names(mraster)[metric]
          
          #--- Character ---#
        } else if (is.character(metric)) {
          if (!metric %in% names(mraster)) {
            stop(glue::glue("'mraster' must have an attribute named {metric}"))
          }
          
        }
      }
    }
    
    message(glue::glue('Implementing optimal allocation of samples based on variability of {metric}'))
    
    #--- merge sraster and mraster metric together ---#
    
    r <- c(sraster,mraster[metric])
    
    vals <- terra::values(r) %>%
      as.data.frame() %>%
      dplyr::select(strata, !!as.name(metric)) %>%
      dplyr::filter(complete.cases(.)) %>%
      dplyr::group_by(strata)
    
    #--- determine number of samples within each strata -- optimal allocation method ---#
    toSample <- vals %>%
      dplyr::summarize(sd = sd(!!as.name(metric)),
                       count = dplyr::n()) %>%
      dplyr::mutate(denom = sum(count*sd)) %>%
      dplyr::rowwise() %>%
      #--- optimal allocation equation ---#
      dplyr::mutate(total = round((nSamp*count*sd)/denom)) %>%
      dplyr::select(strata, total)
    
    
  }
  
  if(allocation == "equal"){
    
    #--- error handling when allocation algorithm is 'unique' ---#
    
    vals <- terra::values(sraster) %>%
      as.data.frame() %>%
      stats::na.omit()
    
    #--- determine total strata ---#
    
    totStrat <- length(unique(vals$strata))
    
    if (!(nSamp %% totStrat == 0)) {
      stop("allocation = 'equal' - nSamp must be divisible by number of strata in 'sraster'.")
    }
    
    
    toSample <- vals %>%
      group_by(strata) %>%
      dplyr::summarize(count = dplyr::n(),
                       total = nSamp / totStrat)
    
    
  }
  
  #--- calculate total samples allocated ---#
  
  tot <- sum(toSample$total)
  
  #--- determine whether there is a difference between 'nSamp' and the number of allocated samples with each stratum ---#
  
  diff <- tot - nSamp
  
  #--- if there is a difference and existing samples are not included ---#
  if(!missing(existing)){
    
    #--- if existing is provided include already sampled plots to achieve the total number ---#
    
    if (!inherits(existing, "data.frame") && !inherits(existing, "sf")) {
      stop("'existing' must be a data.frame or sf object")
    }
    
    if (any(!c("strata") %in% names(existing))) {
      stop("'existing' must have an attribute named 'strata'. Consider using extract_strata().")
    }
    
    #--- convert existing to data frame of strata values ---#
    
    existing <- data.frame(strata = existing$strata)
    
    #--- determine number of samples for each strata ---#
    
    existing <- existing %>%
      dplyr::group_by(strata) %>%
      dplyr::arrange() %>%
      dplyr::summarize(eTotal = dplyr::n())
    
    #--- if the strata for toSample and existing are not identical throw an error ---#
    if (!identical(unique(existing$strata), unique(toSample$strata))) {
      stop("Strata for 'sraster' and 'existing' are not identical. Consider using extract_strata().")
    }
    
    #--- join the 2 df together and subtract the number of existing plots by strata from toSample ---#
    toSample <- toSample %>%
      dplyr::left_join(existing, by = "strata") %>%
      dplyr::mutate(total = total - eTotal) %>%
      dplyr::select(-eTotal) %>%
      as.data.frame()
    
    toSample$need <- existing$eTotal + toSample$total
    
    #--- calculate total samples allocated ---#
    
    tot <- sum(toSample$need)
    
    #--- determine whether there is a difference between 'nSamp' and the number of allocated samples with each stratum ---#
    
    diff <- tot - nSamp
  }
  
  if (diff != 0) {
    
    #--- adjust sample count to force the user defined number ---#
    
    if(force == TRUE){
      
      message(glue::glue('Forcing {nSamp} total samples.'))
      
      #--- if samples need to removed ---#
      
      if (diff > 0){
        
        diffAbs <- abs(diff)
        
        while(diffAbs > 0){
          
          stratAdd <- toSample %>%
            {if (nrow(dplyr::filter(toSample, total == max(total))) > 0) as.data.frame(dplyr::filter(toSample,total == max(total))) else as.data.frame(filter(toSample,total < max(total)))} %>%
            dplyr::sample_n(1) %>%
            dplyr::select(strata) %>%
            dplyr::pull()
          
          toSample <- toSample %>%
            dplyr::mutate(total=replace(total, strata==stratAdd, total[strata==stratAdd] - 1))
          
          diffAbs <- diffAbs - 1
          
        }
        
        #--- if samples need to added ---#
        
      } else if (diff < 0){
        
        diffAbs <- abs(diff)
        
        while(diffAbs > 0){
          
          stratAdd <- toSample %>%
            {if (nrow(dplyr::filter(toSample, total == min(total))) > 0) as.data.frame(dplyr::filter(toSample,total == min(total))) else as.data.frame(filter(toSample,total > min(total)))} %>%
            dplyr::sample_n(1) %>%
            dplyr::select(strata) %>%
            dplyr::pull()
          
          toSample <- toSample %>%
            dplyr::mutate(total=replace(total, strata==stratAdd, total[strata==stratAdd] + 1))
          
          diffAbs <- diffAbs - 1
          
        }
        
      }
      
    } else {
      
      
      message(glue::glue('nSamp of {nSamp} is not perfectly divisible based on strata distribution. nSamp of {tot} will be returned. Use "force = TRUE" to brute force to {nSamp}.'))
    }
  }
  
  toSample
}

