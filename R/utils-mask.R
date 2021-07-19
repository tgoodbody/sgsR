#' Masking
#'
#' Create covariate and sample matrices
#'
#' @inheritParams sample_srs
#' @param raster SpatRaster. Raster to be masked.
#' @family masking
#' @name masking
NULL


#' @export
#' @rdname masking
#' @family masking
#' @keywords internal

mask_access <- function(raster,
                        access,
                        buff_inner,
                        buff_outer) {

  #--- list all buffers to catch NULL values within error handling ---#
  buffers <- list(buff_inner, buff_outer)

  #--- error handling in the presence of 'access' ---#
  if (any(vapply(buffers, is.null, TRUE))) {
    stop("All 'buff_*' paramaters must be provided when 'access' is defined.")
  }

  if (!any(vapply(buffers, is.numeric, FALSE))) {
    stop("All 'buff_*' paramaters must be type numeric")
  }

  message(
    paste0(
      "An access layer has been provided. An internal buffer of ",
      buff_inner,
      " m and an external buffer of ",
      buff_outer,
      " m have been applied"
    )
  )

  #--- convert vectors to spatVector to synergize with terra raster functions---#
  roads <- terra::vect(access)

  #--- make access buffer with user defined values ---#

  buff_in <- terra::buffer(
    x = roads,
    width = buff_inner
  )

  buff_out <- terra::buffer(
    x = roads,
    width = buff_outer
  )

  #--- make difference and aggregate inner and outer buffers to prevent sampling too close to access ---#
  buffer <- terra::aggregate(buff_out - buff_in)

  raster <- terra::mask(raster, mask = buffer)

  #--- make a list with the raster and buffer as outputs ---#

  lout <- list(rast = raster, buff = buffer)

  #--- return the out list ---#

  return(lout)
}
