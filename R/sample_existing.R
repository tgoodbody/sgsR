#' Sample existing
#'
#' @description Sub-sample an existing sample. Four sampling methods are available:
#' \code{clhs}, \code{balanced}, \code{srs} and \code{strat}.
#'
#' @family sample functions
#'
#' @inheritParams sample_systematic
#' @inheritParams extract_strata
#' @inheritParams sample_clhs
#'
#' @param raster SpatRaster. Raster to guide the location of the samples. If \code{type = "clhs"} this raster can also
#' be used to define the population distributions to be used for sampling.
#' @param type Character. A string indicating the type of sampling method to use.
#' Possible values are \code{"clhs"}, \code{"balanced"}, \code{"srs"} and \code{"strat"}.
#' @param ... Additional arguments for the sampling method selected.
#'
#' @return An sf object of samples or a list object if \code{details = TRUE}
#'
#' @note When \code{type = "clhs"} or \code{type = "balanced"} all attributes in \code{existing} will be used for sampling.
#' Remove attributes not indented for sampling' prior to using this algorithm.
#'
#' @author Tristan R.H. Goodbody
#'
#' @export

sample_existing <- function(existing,
                            nSamp,
                            raster = NULL,
                            type = "clhs",
                            access = NULL,
                            buff_inner = NULL,
                            buff_outer = NULL,
                            plot = FALSE,
                            details = FALSE,
                            filename = NULL,
                            overwrite = FALSE,
                            ...) {
  #--- error handling ---#
  
  if(!type %in% c("clhs","balanced","srs","strat")){
    stop("'type' must be one of 'clhs','balanced', 'srs', 'strat'.", call. = FALSE)
  }

  check_existing(
    existing = existing,
    raster = raster,
    nSamp = nSamp,
    plot = plot,
    details = details
  )

  existing <- prepare_existing(
    existing = existing,
    raster = raster,
    access = access,
    buff_inner = buff_inner,
    buff_outer = buff_outer
  )

  #--- sampling ---#

  if (type == "clhs") {
    samples <- sample_existing_clhs(
      existing = existing,
      nSamp = nSamp,
      filename = filename,
      details = details,
      overwrite = overwrite,
      raster = raster,
      ...
    )
  }

  if (type == "balanced") {
    samples <- sample_existing_balanced(
      existing = existing,
      nSamp = nSamp,
      filename = filename,
      overwrite = overwrite,
      ...
    )
  }

  if (type == "srs") {
    samples <- sample_existing_srs(
      existing = existing,
      nSamp = nSamp,
      filename = filename,
      overwrite = overwrite
    )
  }

  if (type == "strat") {
    toSample <- calculate_allocation_existing(
      existing = existing,
      nSamp = nSamp,
      ...
    )

    samples <- sample_existing_strat(
      existing = existing,
      toSample = toSample,
      filename = filename,
      overwrite = overwrite
    )

    if (isTRUE(details)) {
      samples <- list(samples = samples, details = toSample)
    }
  }

  return(samples)
}
