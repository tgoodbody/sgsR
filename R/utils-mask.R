#' Masking
#'
#' Mask covariates and existing samples by access
#'
#' @inheritParams sample_srs
#' @inheritParams extract_strata
#' @param raster SpatRaster. Raster to be masked.
#' @family masking
#' @name masking
#' @return Raster/existing samples masked by provided \code{access} layer and buffers.
NULL


#' @export
#' @rdname masking
#' @family masking
#' @keywords internal

mask_access <- function(raster,
                        access,
                        buff_inner = NULL,
                        buff_outer) {
  #--- error handling in the presence of 'access' ---#
  if (!inherits(access, "sf")) {
    stop("'access' must be an 'sf' object.", call. = FALSE)
  }

  if (!inherits(sf::st_geometry(access), "sfc_MULTILINESTRING") && !inherits(sf::st_geometry(access), "sfc_LINESTRING")) {
    stop("'access' geometry type must be 'LINESTRING' or 'MULTILINESTRING'.", call. = FALSE)
  }

  #--- only specifying an external buffer ---#

  if (is.null(buff_inner)) {
    #--- error handling in the presence of 'access' ---#
    if (is.null(buff_outer)) {
      stop("'buff_outer' must be provided when 'access' is defined.", call. = FALSE)
    }

    if (!is.numeric(buff_outer)) {
      stop("'buff_outer' must be type numeric.", call. = FALSE)
    }

    message(
      paste0("An access layer has been provided. An external buffer of ", buff_outer, " m have been applied.")
    )

    #--- convert vectors to spatVector to synergize with terra raster functions---#
    roads <- terra::vect(access)

    #--- make access buffer with user defined values ---#

    buff_out <- terra::aggregate(
      terra::buffer(
        x = roads,
        width = buff_outer
      )
    )

    raster <- terra::mask(raster, mask = buff_out)

    #--- make a list with the raster and buffer as outputs ---#

    lout <- list(rast = raster, buff = buff_out)

    #--- return the out list ---#

    return(lout)
  } else {
    #--- specifying both internal and external buffers ---#

    #--- list all buffers to catch NULL values within error handling ---#
    buffers <- list(buff_inner, buff_outer)

    #--- error handling in the presence of 'access' ---#

    if (buff_inner > buff_outer) {
      stop("'buff_inner' must be < 'buff_outer'.", call. = FALSE)
    }

    if (any(vapply(buffers, is.null, TRUE))) {
      stop("All 'buff_*' paramaters must be provided when 'access' is defined.", call. = FALSE)
    }

    if (!any(vapply(buffers, is.numeric, FALSE))) {
      stop("All 'buff_*' paramaters must be type numeric.", call. = FALSE)
    }

    message(
      paste0("An access layer has been provided. An internal buffer of ", buff_inner, " m and an external buffer of ", buff_outer, " m have been applied.")
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
}

#' @export
#' @rdname masking
#' @family masking
#' @keywords internal

mask_existing <- function(access,
                          existing,
                          buff_inner = NULL,
                          buff_outer) {
  #--- set global variables ---#

  accessconstraint <- NULL

  #--- error handling in the presence of 'access' ---#
  if (!inherits(access, "sf")) {
    stop("'access' must be an 'sf' object.", call. = FALSE)
  }

  if (!inherits(sf::st_geometry(access), "sfc_MULTILINESTRING") && !inherits(sf::st_geometry(access), "sfc_LINESTRING")) {
    stop("'access' geometry type must be 'LINESTRING' or 'MULTILINESTRING'.", call. = FALSE)
  }

  if (isFALSE(sf::st_crs(access) == sf::st_crs(existing))) {
    stop("'access' and 'existing' do not have the same coordinate system.", call. = FALSE)
  }

  #--- only specifying an external buffer ---#

  if (is.null(buff_inner)) {
    #--- error handling in the presence of 'access' ---#
    if (is.null(buff_outer)) {
      stop("'buff_outer' must be provided when 'access' is defined.", call. = FALSE)
    }

    if (!is.numeric(buff_outer)) {
      stop("'buff_outer' must be type numeric.", call. = FALSE)
    }

    message(
      paste0("An access layer has been provided. An external buffer of ", buff_outer, " m have been applied.")
    )

    #--- buffer access and take all intersecting sample units ---#
    buffer <- sf::st_buffer(x = access, dist = buff_outer) %>%
      sf::st_union(.) %>%
      sf::st_as_sf(.)


    samples <- buffer %>%
      dplyr::mutate(accessconstraint = 1) %>%
      sf::st_join(existing, .) %>%
      dplyr::filter(!is.na(accessconstraint)) %>%
      dplyr::select(-accessconstraint)
  } else {
    #--- specifying both internal and external buffers ---#

    #--- list all buffers to catch NULL values within error handling ---#
    buffers <- list(buff_inner, buff_outer)

    #--- error handling in the presence of 'access' ---#

    if (buff_inner > buff_outer) {
      stop("'buff_inner' must be < 'buff_outer'.", call. = FALSE)
    }

    if (any(vapply(buffers, is.null, TRUE))) {
      stop("All 'buff_*' paramaters must be provided when 'access' is defined.", call. = FALSE)
    }

    if (!any(vapply(buffers, is.numeric, FALSE))) {
      stop("All 'buff_*' paramaters must be type numeric.", call. = FALSE)
    }

    message(
      paste0("An access layer has been provided. An internal buffer of ", buff_inner, " m and an external buffer of ", buff_outer, " m have been applied.")
    )

    #--- convert vectors to spatVector to synergize with terra raster functions---#
    buff_indv <- lapply(X = buffers, function(x) sf::st_buffer(x = access, dist = x) %>% sf::st_union(.))

    buffer <- sf::st_difference(x = buff_indv[[2]], y = buff_indv[[1]]) %>%
      sf::st_as_sf(.)

    samples <- buffer %>%
      dplyr::mutate(accessconstraint = 1) %>%
      sf::st_join(existing, .) %>%
      dplyr::filter(!is.na(accessconstraint)) %>%
      dplyr::select(-accessconstraint)
  }

  message(
    paste0("Masking resulted in an output of ", nrow(samples), " potential sample units.")
  )

  return(list(samples = samples, buffer = buffer))
}
