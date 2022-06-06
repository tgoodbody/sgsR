#' Sample allocation type and count
#'
#' @description Determine how many samples to allocate within strata.
#'
#' @family calculate functions
#'
#' @inheritParams sample_srs
#' @inheritParams sample_strat
#' @param allocation Character. Allocation algorithm to be used. Either \code{prop} (default) for proportional allocation,
#'  \code{optim} for optimal allocation (equal sampling cost), \code{equal} for equal number of samples (defined by \code{nSamp})
#'  for each strata, or \code{"manual"} for user defined strata weights defined using \code{weights}.
#' @param weights Numeric. Only applicable when \code{allocation = "manual"}. Vector of weights where \code{sum(weights) == 1}. Vector length
#' must be equal to the number of unique strata where the first numeric value corresponds to stratum 1, second stratum 2 etc. 
#' @param mraster spatRaster. ALS metric raster. Required when \code{allocation = optim}.
#' @param force Logical. \code{Default = FALSE} - force \code{nSamp} to be exactly the user defined value
#' in cases where \code{nSamp} and \code{sraster} strata count are not equally divisible. Additional samples often need to be allocated or removed
#' based on rounding differences resulting from proportional differences between \code{nSamp} and strata coverages in \code{sraster}.
#' In these instances samples are either added to strata with the lowest number of samples or are removed from strata with the highest number of samples.
#' Has no effect when \code{existing} is provided.
#'
#' @return Returns a data.frame of:
#' \itemize{
#' \item{strata} - Strata ID.
#' \item{total} - Total samples to be allocated based on under representation (positive) or over representation (negative)
#' to be removed at the users discretion.
#' \item{need} - Total required samples per strata. Rounded.
#' }
#'
#' @references
#' Gregoire, T.G., & Valentine, H.T. (2007). Sampling Strategies for Natural Resources and the Environment (1st ed.).
#'  Chapman and Hall/CRC. https://doi.org/10.1201/9780203498880
#'
#' @examples
#' #--- Load strata raster and existing samples---#
#' r <- system.file("extdata", "sraster.tif", package = "sgsR")
#' sr <- terra::rast(r)
#'
#' e <- system.file("extdata", "existing.shp", package = "sgsR")
#' e <- sf::st_read(e)
#'
#' #--- proportional allocation ---#
#' calculate_allocation(
#'   sraster = sr,
#'   nSamp = 50
#' )
#' 
#' #--- equal allocation ---#
#' calculate_allocation(
#'   sraster = sr,
#'   allocation = "equal",
#'   nSamp = 10
#' )
#' 
#' #--- manual allocation ---#
#' #--- define user-defined weights ---#
#' 
#' weights <- c(0.2, 0.2, 0.5, 0.1)
#' 
#' calculate_allocation(
#'   sraster = sr,
#'   allocation = "manual",
#'   weights = weights,
#'   nSamp = 200
#' )
#'
#' @author Tristan R.H. Goodbody
#'
#' @export

calculate_allocation <- function(sraster,
                                 nSamp,
                                 allocation = "prop",
                                 weights = NULL,
                                 mraster = NULL,
                                 existing = NULL,
                                 force = FALSE) {

  #--- Error management ---#
  if (!inherits(sraster, "SpatRaster")) {
    stop("'sraster' must be type SpatRaster", call. = FALSE)
  }
  
  if (any(!c("strata") %in% names(sraster))) {
    stop("'sraster must have a layer named 'strata'", call. = FALSE)
  }
  
  #--- if the sraster has multiple bands subset the band named strata ---#
  if (terra::nlyr(sraster) > 1) {
    stop("Multiple layers detected in sraster. Provide only a single band.", call. = FALSE)
  }

  if (!is.numeric(nSamp)) {
    stop("'nSamp' must be type numeric", call. = FALSE)
  }

  if (!any(allocation == c("prop", "optim", "equal", "manual"))){
    stop(paste0("Unknown allocation type: '", allocation,"' provided. Please use 'prop' (default), 'optim', 'equal', or 'manual'."), call. = FALSE)
  }

  if (!is.logical(force)) {
    stop("'force' must be type logical", call. = FALSE)
  }

  #--- set global vars ---#

  strata <- total <- eTotal <- NULL

  #--- determine crs of input sraster ---#
  crs <- terra::crs(sraster, proj = TRUE)
  
  #--- determine which allocation algorithm to use ---#

  if (allocation != "equal") {
    
    #--- proportional allocation ---#
    
    if (allocation == "prop") {
      
      #--- error handling when allocation algorithm is 'prop' ---#

      if (!is.null(mraster)) {
        message("'mraster' was specified but 'allocation = prop' - did you mean to use 'allocation = optim'?")
      }
      
      if (!is.null(weights)) {
        message("'weights' was specified but 'allocation = prop' - did you mean to use 'allocation = manual'?")
      }

      toSample <- allocation_prop(sraster = sraster, nSamp = nSamp)
    }

    #--- optimal allocation ---#

    if (allocation == "optim") {

      #--- error handling when allocation algorithm is 'optim' ---#

      if (!is.null(weights)) {
        message("'weights' was specified but 'allocation = optim' - did you mean to use 'allocation = manual'?")
      }

      toSample <- allocation_optim(sraster = sraster, mraster = mraster, nSamp = nSamp)
      
    }
    
    #--- manual allocation ---#
    if (allocation == "manual"){
      
      #--- error handling when allocation algorithm is 'manual' ---#
      
      if (!is.null(mraster)) {
        message("'allocation = manual' - ignoring 'mraster'")
      }
      
      toSample <- allocation_manual(sraster = sraster, nSamp = nSamp, weights = weights)
    }

    #--- calculate total samples allocated ---#

    tot <- sum(toSample$total)
    
  } else {
    
    #--- equal allocation ---#
    
    if (!is.null(weights)) {
      message("'weights' was specified but 'allocation = equal' - did you mean to use 'allocation = manual'?")
    }
    
    if (!is.null(mraster)) {
      message("'mraster' was specified but 'allocation = equal' - did you mean to use 'allocation = optim'?")
    }

    toSample <- allocation_equal(sraster = sraster, nSamp = nSamp)

    tot <- unique(toSample$total)
  }

  #--- determine whether there is a difference between 'nSamp' and the number of allocated samples with each stratum ---#

  diff <- tot - nSamp

  #--- if there is a difference and existing samples are not included ---#
  if (!missing(existing)) {

    #--- if existing is provided include already sampled plots to achieve the total number ---#

    if (!inherits(existing, "data.frame") && !inherits(existing, "sf")) {
      stop("'existing' must be a data.frame or sf object", call. = FALSE)
    }

    if (any(!c("strata") %in% names(existing))) {
      stop("'existing' must have an attribute named 'strata'. Consider using extract_strata().", call. = FALSE)
    }

    #--- convert existing to data frame of strata values ---#

    existing <- data.frame(strata = existing$strata)

    #--- determine number of samples for each strata ---#

    existing <- existing %>%
      dplyr::group_by(strata) %>%
      dplyr::arrange() %>%
      dplyr::summarize(eTotal = dplyr::n())
    
    #--- check if samples fall in areas where stratum values are NA ---#
    
    if(any(!complete.cases(existing$strata))){
      
      nNA <- existing %>%
        dplyr::filter(!complete.cases(strata)) %>%
        dplyr::pull(eTotal)
      
      message(paste0(nNA," samples in `existing` are located where strata values are NA."))

      existing <- existing %>%
        stats::na.omit()
      
    }

    #--- if the strata for toSample and existing are not identical throw an error ---#
    if (!all.equal(unique(existing$strata), unique(toSample$strata))) {
      stop("Strata for 'sraster' and 'existing' are not identical. Consider using extract_strata().", call. = FALSE)
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

  if (allocation != "equal") {
    if (diff != 0) {

      #--- adjust sample count to force the user defined number ---#

      if (force == TRUE) {
        message(paste0("Forcing ", nSamp, " total samples."))

        #--- if samples need to be removed ---#

        if (diff > 0) {
          diffAbs <- abs(diff)

          while (diffAbs > 0) {
            stratAdd <- toSample %>%
              {
                if (nrow(dplyr::filter(toSample, total == max(total))) > 0) as.data.frame(dplyr::filter(toSample, total == max(total))) else as.data.frame(dplyr::filter(toSample, total < max(total)))
              } %>%
              dplyr::sample_n(1) %>%
              dplyr::select(strata) %>%
              dplyr::pull()

            toSample <- toSample %>%
              dplyr::mutate(total = replace(total, strata == stratAdd, total[strata == stratAdd] - 1))

            diffAbs <- diffAbs - 1
          }

          #--- if samples need to be added ---#
        } else if (diff < 0) {
          diffAbs <- abs(diff)

          while (diffAbs > 0) {
            stratAdd <- toSample %>%
              {
                if (nrow(dplyr::filter(toSample, total == min(total))) > 0) as.data.frame(dplyr::filter(toSample, total == min(total))) else as.data.frame(dplyr::filter(toSample, total > min(total)))
              } %>%
              dplyr::sample_n(1) %>%
              dplyr::select(strata) %>%
              dplyr::pull()

            toSample <- toSample %>%
              dplyr::mutate(total = replace(total, 
                                            strata == stratAdd, 
                                            total[strata == stratAdd] + 1)
                            )

            diffAbs <- diffAbs - 1
          }
        }
      } else {
        message(paste0("nSamp of ",nSamp," is not perfectly divisible based on strata distribution. nSamp of ", tot, " will be returned. Use 'force = TRUE' to brute force to ", nSamp,"."))
      }
    }
  } else {
    if (force == TRUE) {
      message("force == TRUE has no effect when allocation == 'equal'. Ignorning")
    }
  }

  toSample
}
