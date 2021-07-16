#' Extract
#'
#' Extract raster values to samples
#'
#' @family extract
#' @inheritParams strat_breaks
#' @inheritParams strat_kmeans
#' @param existing sf. Samples resulting from sample_* functions.
#' @param data.frame Logical. If true outputs as data.frame
#' @name extract
#' @importFrom magrittr %>%
#' @importFrom methods is
NULL

#' @export
#' @rdname extract
#' @family extract
#' @param sraster spatRaster. Stratification raster.
#' @keywords internal
#' @return An sf or data.frame object of samples with strata attributes
#' #--- function for Extract and buffering access ---#


extract_strata <- function(sraster,
                           existing,
                           data.frame = FALSE) {

  #--- Set global vars ---#

  x <- y <- X <- Y <- strata <- NULL

  #--- Error management ---#

  if (!inherits(sraster, "SpatRaster")) {
    stop("'sraster' must be type SpatRaster", call. = FALSE)
  }

  if (any(!c("strata") %in% names(sraster))) {
    stop("'sraster' must have a layer named 'strata'")
  }

  if (!inherits(existing, "sf") && inherits(sf::st_geometry(existing), "sfc_POINT")) {
    stop("'existing' must be an 'sf' object of type 'sfc_POINT' geometry")
  }

  if (!is(existing, "data.frame")) {
    stop("existing must be a data.frame")
  }

  #--- if the existing plots are an sf object extract coordinates ---#

  if (is(existing, "sf")) {

    #--- Convert to spatVector to enable extraction of strata values ---#

    existing <- sf::st_coordinates(existing)
  } else {
    if (any(!c("X", "Y") %in% colnames(existing))) {

      #--- if coordinate column names are lowercase change them to uppercase to match requirements ---#

      if (any(c("x", "y") %in% colnames(existing))) {
        existing <- existing %>%
          dplyr::rename(
            X = x,
            Y = y
          )

        message("Column coordinate names are lowercase - converting to uppercase")
      } else {

        #--- if no x/y columns are present stop ---#

        stop("'existing' must have columns named 'X' and 'Y'")
      }
    }
  }

  #--- extract values from the sraster dataset ---#

  strata_vals <- terra::extract(sraster, existing)

  #--- bind values and coordinates ---#

  existing_strata <- cbind(existing, strata_vals)

  #--- select only coordinate and strata values ---#

  existing_strata <- existing_strata %>%
    dplyr::select(X, Y, strata)

  #--- output either data.frame or sf object ---#

  if (isTRUE(data.frame)) {

    #--- return data.frame ---#
    return(existing_strata)
  } else {

    #--- convert coordinates to a sf object ---#

    samples <- existing_strata %>%
      as.data.frame() %>%
      sf::st_as_sf(., coords = c("X", "Y"))

    #--- assign sraster crs to spatial points object ---#

    sf::st_crs(samples) <- terra::crs(sraster)

    #--- return sf object ---#
    return(samples)
  }
}


#' @export
#' @rdname extract
#' @family extract
#' @description Extract metric raster attributes to existing
#' @return An sf or data.frame object of samples with associated raster cell attributes
#' #--- Extract raster metrics to existing sample ---#

extract_metrics <- function(mraster,
                            existing,
                            data.frame = FALSE) {

  #--- Set global vars ---#

  x <- y <- X <- Y <- strata <- geometry <- NULL

  #--- Error management ---#

  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  }

  if (!inherits(existing, "sf")) {
    stop("'existing' must be an 'sf' object", call. = FALSE)
  }

  #--- Extract coordinates from existing ---#

  xy <- sf::st_coordinates(existing)

  vals <- terra::extract(mraster, xy)

  #--- extract other attributes from sampling and remove geometry attribute ---#

  samp_mets <- as.data.frame(existing)

  samp_mets <- dplyr::select(samp_mets, -geometry)

  #--- bind values and coordinates ---#
  samples <- cbind(xy, samp_mets, vals)

  if (isTRUE(data.frame)) {

    #--- return data.frame ---#
    return(samples)
  } else {

    #--- convert coordinates to a sf object ---#

    samples <- samples %>%
      as.data.frame() %>%
      sf::st_as_sf(., coords = c("X", "Y"))

    #--- assign mraster crs to spatial points object ---#

    sf::st_crs(samples) <- terra::crs(mraster)

    #--- return sf object ---#
    return(samples)
  }
}
