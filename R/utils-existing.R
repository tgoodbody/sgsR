#' Check existing sample data against requirements
#'
#' This function checks whether the existing sample data meets certain requirements for use in downstream analyses.
#'
#' @inheritParams sample_existing
#'
#' @return If requirements are met, the function returns the prepared existing sample data.
#' Otherwise, it raises a stop error with a relevant message.
#'
#' @keywords internal
#' @noRd
check_existing <- function(existing,
                           raster,
                           nSamp,
                           plot = FALSE,
                           details = NULL) {
  x <- y <- NULL

  #--- error handling ---#

  if (!inherits(existing, "data.frame") && !inherits(existing, "sf")) {
    stop("'existing' must be a data.frame or sf object.", call. = FALSE)
  }

  if (!is.null(raster) & !inherits(raster, "SpatRaster")) {
    stop("'raster' must be type SpatRaster.", call. = FALSE)
  }

  if (!is.numeric(nSamp)) {
    stop("'nSamp' must be type numeric of type integer or double.", call. = FALSE)
  }

  if (nSamp >= nrow(existing)) {
    stop("'nSamp' must be less than the total number of 'existing' samples.", call. = FALSE)
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical.", call. = FALSE)
  }

  if (!is.null(details)) {
    if (!is.logical(details)) {
      stop("'details' must be type logical.", call. = FALSE)
    }
  }

  #--- Prepare existing sample data ---#
  if (!inherits(existing, "sf")) {
    if (any(!c("X", "Y") %in% colnames(existing))) {
      #--- if coordinate column names are lowercase change them to uppercase to match requirements ---#

      if (any(c("x", "y") %in% colnames(existing))) {
        existing <- existing %>%
          dplyr::rename(
            X = x,
            Y = y
          )

        message("'existing' column coordinate names are lowercase - converting to uppercase.")
      } else {
        #--- if no x/y columns are present stop ---#

        stop("'existing' must have columns named 'X' and 'Y'.", call. = FALSE)
      }
    }
  }
}
#' Prepare existing sample data
#'
#' This function prepares the existing sample data by ensuring that it meets the necessary requirements for downstream analysis.
#' If the 'existing' object is not of class 'sf', this function checks that the data contains columns named "X" and "Y", and
#' converts lowercase "x" and "y" column names to uppercase if necessary. If the 'raster' object is supplied, this function
#' checks if 'existing' contains attributes with the same names as 'raster'. If it does not, this function extracts metrics
#' at existing sample locations using the 'raster' object. If 'access' is not null, the function masks the existing sample
#' locations with a buffered access area.
#'
#' @inheritParams sample_existing
#'
#' @keywords internal
#' @noRd
prepare_existing <- function(existing,
                             raster = NULL,
                             access = NULL,
                             buff_inner = NULL,
                             buff_outer = NULL) {
  x <- y <- NULL

  #--- Prepare existing sample data ---#
  if (!inherits(existing, "sf")) {
    if (any(!c("X", "Y") %in% colnames(existing))) {
      #--- if coordinate column names are lowercase change them to uppercase to match requirements ---#

      if (any(c("x", "y") %in% colnames(existing))) {
        existing <- existing %>%
          dplyr::rename(
            X = x,
            Y = y
          )

        message("'existing' column coordinate names are lowercase - converting to uppercase.")
      } else {
        #--- if no x/y columns are present stop ---#

        stop("'existing' must have columns named 'X' and 'Y'.", call. = FALSE)
      }
    }

    #--- if raster is supplied
    if (!is.null(raster)) {
      crs <- terra::crs(raster)

      #--- check to see if 'existing' does not contain attributes with the same names as 'raster' ---#
      if (!all(names(raster) %in% names(existing))) {
        message("'existing' does not contain attributes with the same names as 'raster'. Extracting metrics.")

        existing <- existing %>%
          sf::st_as_sf(., coords = c("X", "Y"), crs = crs)

        #--- access constraint ---#

        if (!is.null(access)) {
          masked <- mask_existing(access = access, existing = existing, buff_inner = buff_inner, buff_outer = buff_outer)

          existing <- masked$samples
        }

        #--- extract covariates at existing sample locations ---#
        existing <- extract_metrics(mraster = raster, existing = existing) %>%
          na.omit()
      }
    }

    existing <- existing
  } else {
    #---  bring existing crs forward ---#
    crs <- sf::st_crs(existing)

    #--- check to see if 'existing' does not contain attributes with the same names as 'raster' ---#
    if (!all(names(raster) %in% names(existing))) {
      message("'existing' does not contain attributes with the same names as 'raster'. Extracting metrics.")

      #--- access constraint ---#

      if (!is.null(access)) {
        masked <- mask_existing(access = access, existing = existing, buff_inner = buff_inner, buff_outer = buff_outer)

        existing <- masked$samples
      }

      #--- extract covariates at existing sample locations ---#
      existing <- extract_metrics(mraster = raster, existing = existing) %>%
        na.omit()
    }
  }

  return(existing)
}

#' Get existing and XY coordinates
#'
#' This function extracts the 'X' and 'Y' coordinates from an 'sf' object and returns them in a data.frame format.
#'
#' @inheritParams sample_existing
#'
#' @return \code{existing} with respective 'X' and 'Y' columns representing the coordinates of the input \code{sf} object.
#'
#' @keywords internal
#' @noRd
coords_existing <- function(existing) {
  xy <- existing %>%
    sf::st_coordinates(.)

  existing <- existing %>%
    sf::st_drop_geometry(.) %>%
    cbind(xy, .)

  existing
}
