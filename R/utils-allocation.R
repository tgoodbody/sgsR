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

allocation_prop <- function(sraster,
                            nSamp){
  
  #--- define global vars ---#
  
  strata <- count <- freq <- total <- NULL
  
  message("Implementing porportional allocation of samples")
  
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

allocation_optim <- function(sraster,
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
  
  message(paste0("Implementing optimal allocation of samples based on variability of '", nm,"'"))
  
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

allocation_manual <- function(sraster,
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
  
  message("Implementing allocation of samples based on user-defined weights")

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

allocation_equal <- function(sraster,
                             nSamp){
  
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