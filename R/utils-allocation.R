#' Allocating
#'
#' Allocation algorithms
#'
#' @inheritParams calculate_allocation
#' @family allocating
#' @name allocating
NULL


#' @export
#' @rdname allocating
#' @family allocation
#' @keywords internal

allocate_prop <- function(sraster,
                          nSamp){
  
  #--- define global vars ---#
  
  strata <- count <- freq <- total <- NULL
  
  message("Implementing proportional allocation of samples.")
  
  #--- generate vals data.frame ---#
  
  vals <- terra::values(sraster) %>%
    as.data.frame() %>%
    stats::na.omit()

  #--- determine number of samples within each strata ---#
  toSample <- vals %>%
    dplyr::group_by(strata) %>%
    dplyr::summarize(count = dplyr::n()) %>%
    dplyr::mutate(freq = count / sum(count),
                  total = freq * nSamp) %>%
    dplyr::mutate(total = replace(total, total < 1, 1)) %>%
    dplyr::mutate(total = round(total, digits = 0)) %>%
    dplyr::select(strata, total) %>%
    as.data.frame()

  toSample

}

#' @export
#' @rdname allocating
#' @family allocation
#' @keywords internal

allocate_optim <- function(sraster,
                             mraster,
                             nSamp){
  
  #--- define global vars ---#
  
  v_sd <- strata <- count <- total <- denom <- NULL
  
  #--- error handling when allocation algorithm is 'optim' ---#
  
  if (is.null(mraster)) {
    stop("'mraster' must be supplied if 'allocation = optim'.", call. = FALSE)
  }
  
  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  }
  
  #--- if there is only 1 band in mraster use it as default ---#
  
  if (terra::nlyr(mraster) == 1) {
    rastermetric <- mraster
    nm <- names(rastermetric)
  } else {
    stop("Multiple layers detected in 'mraster'. Please define a singular band for allocation.", call. = FALSE)
  }
  
  message(paste0("Implementing optimal allocation of samples based on variability of '", nm,"'."))
  
  #--- merge sraster and mraster together ---#
  
  r <- c(sraster, rastermetric)
  
  vals <- terra::values(r) %>%
    as.data.frame() %>%
    dplyr::select(strata, !!nm) %>%
    stats::na.omit() %>%
    dplyr::group_by(strata)
  
  #--- determine number of samples within each strata -- optimal allocation method ---#
  toSample <- vals %>%
    dplyr::summarize(
      v_sd = sd(eval(as.name(nm))),
      count = dplyr::n()) %>%
    dplyr::mutate(denom = sum(count * v_sd)) %>%
    dplyr::rowwise() %>%
    #--- optimal allocation (equal sampling cost) equation. See Gregoire & Valentine (2007) Section 5.4.4 ---#
    dplyr::mutate(total = round(nSamp * ((count * v_sd) / denom)), digits = 0) %>%
    dplyr::select(strata, total)

  toSample
  
}

#' @export
#' @rdname allocating
#' @family allocation
#' @keywords internal

allocate_manual <- function(sraster,
                              nSamp,
                              weights){
  
  #--- define global vars ---#
  
  strata <- total <- NULL
  
  #--- error handling when allocation algorithm is 'manual' ---#
  
  if(is.null(weights)){
    stop("'weights' must be defined if 'allocation = manual'.", call. = FALSE)
  }
  
  if(!is.numeric(weights)){
    stop("'weights' must be a numeric vector.", call. = FALSE)
  }
  
  if(sum(weights) != 1){
    stop("'weights' must add up to 1.", call. = FALSE)
  }
  
  message("Implementing allocation of samples based on user-defined weights.")

  #--- generate vals data.frame ---#
  
  vals <- terra::values(sraster) %>%
    as.data.frame() %>%
    stats::na.omit()

  if(length(weights) != length(unique(vals$strata))){
    stop("'weights' must be the same length as the number of strata in 'sraster'.", call. = FALSE)
  }
  
  #--- determine number of samples within each strata ---#
  toSample <- vals %>%
    dplyr::group_by(strata) %>%
    dplyr::summarize(count = dplyr::n()) %>%
    dplyr::mutate(weights = weights,
                  total = nSamp * weights) %>%
    dplyr::mutate(total = replace(total, total < 1, 1)) %>%
    dplyr::mutate(total = round(total, digits = 0)) %>%
    dplyr::select(strata, total) %>%
    as.data.frame()

  toSample
  
}

#' @export
#' @rdname allocating
#' @family allocation
#' @keywords internal

allocate_equal <- function(sraster,
                             nSamp){
  
  message("Implementing equal allocation of samples.")
  
  #--- define global vars ---#
  
  strata <- NULL
  
  #--- generate vals data.frame ---#
  
  vals <- terra::values(sraster) %>%
    as.data.frame() %>%
    stats::na.omit()

  #--- assign nSamp to each strata ---#
  
  toSample <- vals %>%
    dplyr::group_by(strata) %>%
    dplyr::summarize(total = nSamp)

  toSample
  
}

#' @export
#' @rdname allocating
#' @family allocation
#' @keywords internal

allocate_existing <- function(toSample,
                                existing){
  
  strata <- total <- eTotal <- NULL
  
  #--- if existing is provided include already sampled plots to achieve the total number ---#
  
  if (!inherits(existing, "data.frame") && !inherits(existing, "sf")) {
    stop("'existing' must be a data.frame or sf object.", call. = FALSE)
  }
  
  if (any(!c("strata") %in% names(existing))) {
    stop("'existing' must have an attribute named 'strata'. Consider using extract_strata().", call. = FALSE)
  }
  
  #--- convert existing to data frame of strata values ---#
  
  existing <- data.frame(strata = existing$strata)
  
  #--- determine number of samples for each strata ---#
  
  existing <- existing %>%
    dplyr::group_by(strata) %>%
    dplyr::tally(name = "eTotal")
  
  #--- check if samples fall in areas where stratum values are NA ---#
  
  if(any(!complete.cases(existing$strata))){
    
    nNA <- existing %>%
      dplyr::filter(!complete.cases(strata)) %>%
      dplyr::pull(eTotal)
    
    message(paste0(nNA," samples in `existing` are located where strata values are NA. Expect ",nNA," additional samples in output."))

    existing <- existing %>%
      stats::na.omit()
  }

  #--- if the unique(existing$strata) %in% unique(toSample$strata) for toSample and existing are not identical throw an error ---#
  if (!any(unique(existing$strata) %in% unique(toSample$strata))) {
      stop("'existing' does not contain matching strata to those in `sraster`. Check strata in both data sets & consider using extract_strata().", call. = FALSE)
  }
  
  #--- join the 2 df together and subtract the number of existing plots by strata from toSample ---#
  toSample <- toSample %>%
    dplyr::left_join(existing, by = "strata") %>%
    replace(is.na(.), 0) %>%
    dplyr::mutate(total = total - eTotal,
                  need = eTotal + total) %>%
    dplyr::select(-eTotal)

  toSample
  
}

#' @export
#' @rdname allocating
#' @family allocation
#' @keywords internal

allocate_force <- function(toSample,
                           nSamp,
                           diff){
  
  total <- strata <- NULL
  
  #--- Force the removal of samples to meet 'nSamp' ---#

  if (diff > 0) {
    diffAbs <- abs(diff)
    
    while (diffAbs > 0) {
      stratAdd <- toSample %>%
        {
          if (nrow(dplyr::filter(toSample, total == max(total))) > 0) as.data.frame(dplyr::filter(toSample, total == max(total))) 
          else as.data.frame(dplyr::filter(toSample, total < max(total)))
        } %>%
        dplyr::sample_n(1) %>%
        dplyr::select(strata) %>%
        dplyr::pull()
      
      toSample <- toSample %>%
        dplyr::mutate(total = replace(total, strata == stratAdd, total[strata == stratAdd] - 1))
      
      diffAbs <- diffAbs - 1
    }
    
    #--- Force the addition of samples to meet 'nSamp' ---#
  } else if (diff < 0) {
    diffAbs <- abs(diff)
    
    while (diffAbs > 0) {
      stratAdd <- toSample %>%
        {
          if (nrow(dplyr::filter(toSample, total == min(total))) > 0) as.data.frame(dplyr::filter(toSample, total == min(total))) 
          else as.data.frame(dplyr::filter(toSample, total > min(total)))
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
  
  toSample
}