#' Vectorization helpers
#'
#' @inheritParams strat_breaks
#' @inheritParams strat_quantiles
#' @family vectorize
#' @name vectorize
NULL

#' Breaks vectorize
#' @family vectorize
#' @rdname vectorize
#' @keywords internal

#--- strat_breaks() vectorization function ---#
calculate_breaks <- function(mraster, breaks){
  
  from <- NULL
  
  breaks <- data.frame(from = c(-Inf, breaks, Inf)) %>%
    dplyr::mutate(
      to = dplyr::lead(from),
      becomes = seq(1:length(from))
    ) %>%
    stats::na.omit() %>%
    as.matrix()
  
  rcl <- terra::classify(mraster, breaks, others = NA)
  names(rcl) <- "strata"
  
  return(rcl)
}

#' Quantile breaks
#' @family vectorize
#' @rdname vectorize
#' @keywords internal

calculate_quantile_breaks <- function(mraster, nStrata){
  
  if (!is.numeric(nStrata)) {
    stop("'nStrata' must be type numeric.", call. = FALSE)
  }
  
  if(any(nStrata < 0)){
    
    stop("`nStrata` must be either a positive integer representing the number of desired classes or a numeric vector of probabilities with values in [0,1]", call. = FALSE)
    
  }
  
  if(length(nStrata) > 1 && any(nStrata > 1)){
    
    stop("`nStrata` must be either a positive integer representing the number of desired classes or a numeric vector of probabilities with values in [0,1]", call. = FALSE)
    
  }
  
  #--- check if nStrata is an integer ---#
  if(length(nStrata) < 2){
    
    Qbreaks <- as.vector(as.numeric(terra::global(mraster, quantile_breaks_integer, nStrata = nStrata)))
    
  } else {
    
    Qbreaks <- c(-Inf, as.vector(as.numeric(terra::global(mraster, quantile_breaks, nStrata = nStrata))), Inf)
    
  } 
  
  olead <- dplyr::lead(Qbreaks) %>%
    tidyr::replace_na(Inf)
  
  ftb <- as.matrix(data.frame(from = Qbreaks, to = olead, becomes = 1:length(Qbreaks)))
  
  rout <- terra::classify(mraster, ftb, others = NA, right = NA)
  
  names(rout) <- "strata"
  
  out <- list(raster = rout, breaks = data.frame(names = rep(names(mraster),length(Qbreaks)), val = Qbreaks))
  
  return(out)
}

#' Quantile vectorize by integer helper
#' @family vectorize
#' @rdname vectorize
#' @keywords internal

#--- strat_quantiles vectorization function ---#

quantile_breaks_integer <- function(mraster, nStrata){
  terra::quantile(mraster, probs = seq(0, 1, length.out = (nStrata + 1)), na.rm=TRUE)
}

#' Quantile vectorize by probabilities vector helper
#' @family vectorize
#' @rdname vectorize
#' @keywords internal

quantile_breaks <- function(mraster, nStrata){
  terra::quantile(mraster, probs = nStrata, na.rm=TRUE)
}