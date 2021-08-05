#' Determine required samples in strata
#'
#' @details Count the required number of samples for each raster strata
#'
#' @family calculate functions
#'
#' @inheritParams sample_srs
#' @inheritParams sample_strat
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
#' calculate_reqSamples(sraster = sr, 
#'                      nSamp = 200)
#' 
#' e.sr <- extract_strata(sraster = sr, 
#'                        existing = e)
#' 
#' calculate_reqSamples(sraster = sr, 
#'                      nSamp = 200, 
#'                      existing = e.sr)
#'                      
#' @author Tristan R.H. Goodbody
#'
#' @export


calculate_reqSamples <- function(sraster,
                                 nSamp,
                                 existing = NULL,
                                 force = FALSE) {

  #--- Error management ---#
  if (!inherits(sraster, "SpatRaster")) {
    stop("'sraster' must be type SpatRaster", call. = FALSE)
  }
  
  if (!is.numeric(nSamp)) {
    stop("'nSamp' must be type numeric")
  }
  
  if (!is.logical(force)) {
    stop("'force' must be type logical")
  }
  
  if(isTRUE(force)){
    if(!missing(existing)){
    message("'force = TRUE' and 'existing' is provided. Ignoring 'force'.")
    }
  }
  
  #--- set global vars ---#

  strata <- count <- freq <- total <- eTotal <- NULL

  #--- determine crs of input sraster ---#
  crs <- terra::crs(sraster, proj=TRUE)

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

  #--- determine whether there is a difference between 'n' and the number of samples provided based on count ---#

  diff <- sum(toSample$total) - nSamp

  #--- if there is a difference remove samples from classes with the most samples ---#
  if(missing(existing)){
    
    if (diff != 0) {
    
      #--- adjust sample count to force the user defined number ---#
      
      if(force == TRUE){
        
        message(paste0("Forcing ", nSamp, " total samples."))
        
        #--- if samples need to removed ---#
      
        if (diff > 0){
          
          diffAbs <- abs(diff)
          
          while(diffAbs > 0){
            
            stratAdd <- toSample %>%
              {if (nrow(dplyr::filter(toSample, total == max(total))) > 0) as.data.frame(dplyr::filter(toSample,total == max(total))) else as.data.frame(filter(toSample,total < max(total)))} %>%
              dplyr::sample_n(1) %>%
              dplyr::select(strata) %>%
              pull()
            
            toSample <- toSample %>%
              mutate(total=replace(total, strata==stratAdd, total[strata==stratAdd] - 1))
            
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
              pull()
            
            toSample <- toSample %>%
              mutate(total=replace(total, strata==stratAdd, total[strata==stratAdd] + 1))
            
            diffAbs <- diffAbs - 1
              
          }
          
        }
        
      } else {
        
        
        message(paste0("nSamp of ", nSamp, 
                       " is not equally divisible into ", max(toSample$strata), 
                       " strata. nSamp of ", sum(toSample$total), "  will be returned.",
                       " Use 'force = TRUE' to brute force to ", nSamp, "."))
      }
    }
    
  } else {
    
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
  }

  toSample
}
