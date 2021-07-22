#' Determine required samples in strata
#'
#' @details Count the required number of samples for each raster strata
#'
#' @family calculate functions
#'
#' @inheritParams sample_srs
#' @inheritParams sample_strat
#'
#'
#' @return data.frame of strata and associates samples
#' 
#' @examples 
#' \dontrun{
#' #--- Load strata raster and existing samples---#
#' r <- system.file("extdata","kmeans.tif", package = "sgsR")
#' sr <- terra::rast(r)
#' 
#' e <- system.file("extdata","existing.shp", package = "sgsR")
#' e <- sf::st_read(e)
#'
#' #--- perform grid sampling ---#
#' calculate_reqSamples(sraster = sr, nSamp = 200)
#' 
#' calculate_reqSamples(sraster = sr, nSamp = 200, existing = e)
#' }
#' @author Tristan R.H. Goodbody 
#' 
#' @export


calculate_reqSamples <- function(sraster,
                                 nSamp,
                                 existing = NULL) {

  #--- set global vars ---#

  strata <- count <- freq <- total <- eTotal <- NULL

  #--- determine crs of input sraster ---#
  crs <- terra::crs(sraster)

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

  if (diff != 0) {
    if (diff > 0) {

      #--- determine the largest total samples size among strata ---#

      maxTotal <- max(toSample$total)

      #--- subtract 'diff' from largest sample size ---#

      toSample <- toSample %>%
        dplyr::mutate(total = replace(
          total,
          total == maxTotal,
          maxTotal - abs(diff)
        ))
    } else {

      #--- determine the largest total samples size among strata ---#

      minTotal <- min(toSample$total)

      #--- add 'diff' to smallest sample size ---#

      toSample <- toSample %>%
        dplyr::mutate(total = replace(
          total,
          total == minTotal,
          minTotal + abs(diff)
        ))
    }
  }

  #--- if existing is provided include already sampled plots to achieve the total number ---#

  if (!missing(existing)) {

    #--- convert existing to data frame of strata values ---#

    existing <- data.frame(strata = existing$strata)

    #--- determine number of samples for each strata ---#

    existing <- existing %>%
      dplyr::group_by(strata) %>%
      dplyr::arrange() %>%
      dplyr::summarize(eTotal = dplyr::n())

    #--- if the strata for toSample and existing are not identical throw an error ---#
    if (!identical(unique(existing$strata), unique(toSample$strata))) {
      stop("Strata for 'sraster' and 'existing' are not identical")
    }

    #--- join the 2 df together and subtract the number of existing plots by strata from toSample ---#
    toSample <- toSample %>%
      dplyr::left_join(existing, by = "strata") %>%
      dplyr::mutate(total = total - eTotal) %>%
      dplyr::select(-eTotal) %>%
      as.data.frame()
  }

  toSample
}
