#' Extract metrics
#'
#' @description Extract metric values to existing samples
#'
#' @family extract functions
#'
#' @inheritParams sample_systematic
#' @inheritParams extract_strata
#'
#' @param mraster spatRaster. Metrics Raster.
#' @param data.frame Logical. Output as data.frame if \code{TRUE}
#'
#' @return An sf or data.frame object of samples with metrics attributes.
#'
#' @note
#'
#' If \code{data.frame = TRUE} output will be written using \code{\link[utils]{write.table}}
#'
#' @examples
#' #--- Load mraster ---#
#' r <- system.file("extdata", "mraster.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' #' #--- load existing samples ---#
#' e <- system.file("extdata", "existing.shp", package = "sgsR")
#' e <- sf::st_read(e)
#'
#' extract_metrics(
#'   mraster = mr,
#'   existing = e
#' )
#'
#' @author Tristan R.H. Goodbody
#'
#' @export

extract_metrics <- function(mraster,
                            existing,
                            quiet = FALSE,
                            data.frame = FALSE,
                            filename = NULL,
                            overwrite = FALSE) {
  #--- Set global vars ---#

  x <- y <- X <- Y <- strata <- ID <- geometry <- NULL

  #--- Error management ---#

  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster.", call. = FALSE)
  }

  if (!inherits(existing, "data.frame") && !inherits(existing, "sf")) {
    stop("'existing' must be a data.frame or sf object.", call. = FALSE)
  }

  if (!is.logical(quiet)) {
    stop("'quiet' must be type logical.", call. = FALSE)
  }

  #--- if the existing plots are an sf object extract coordinates ---#

  if (is(existing, "sf")) {
    if (!inherits(sf::st_geometry(existing), "sfc_POINT")) {
      stop("'existing' must be an 'sf' object of type 'sfc_POINT' geometry.", call. = FALSE)
    }

    #--- to preserve input CRS ---#

    crs <- sf::st_crs(existing)

    #--- Extract xy coordinates to enable extraction of strata values ---#

    xy <- sf::st_coordinates(existing)

    existing <- existing %>%
      sf::st_drop_geometry(.) %>%
      dplyr::select(-X, -Y)
  } else {
    #--- To use raster CRS when existing is a data.frame ---#

    crs <- terra::crs(mraster)

    if (any(!c("X", "Y") %in% colnames(existing))) {
      #--- if coordinate column names are lowercase change them to uppercase to match requirements ---#

      if (any(c("x", "y") %in% colnames(existing))) {
        xy <- existing %>%
          dplyr::rename(
            X = x,
            Y = y
          )

        existing <- existing %>%
          dplyr::select(-x, -y)

        message("Column coordinate names are lowercase - converting to uppercase.")
      } else {
        #--- if no x/y columns are present stop ---#

        stop("'existing' must have columns named 'X' and 'Y'.")
      }
    } else {
      xy <- existing %>%
        dplyr::select(X, Y)


      existing <- existing %>%
        dplyr::select(-X, -Y)
    }
  }

  vals <- terra::extract(mraster, xy)

  #--- check that extractions has produced some values -- if not tell the user ---#

  if (all(!complete.cases(vals))) {
    stop("'existing' only extracts NA values. Ensure that 'existing' overlaps with 'mraster'.", call. = FALSE)
  }

  #--- if existing samples are co-located with NA values ---#
  if (any(!complete.cases(vals))) {
    nNA <- vals %>%
      dplyr::filter(!complete.cases(.)) %>%
      dplyr::tally() %>%
      dplyr::pull()

    if (isFALSE(quiet)) {
      message(paste0(nNA, " samples are located where metric values are NA."))
    }
  }


  if (isTRUE(data.frame)) {
    samples <- cbind(xy, vals, existing)

    if ("geometry" %in% names(samples)) {
      samples <- samples %>%
        dplyr::select(-geometry)
    }

    #--- write outputs if desired ---#
    write_samples_df(samples = samples, filename = filename, overwrite = overwrite)

    #--- return data.frame ---#
    return(samples)
  } else {
    #--- convert coordinates to a sf object ---#

    samples <- cbind(xy, vals, existing) %>%
      sf::st_as_sf(., coords = c("X", "Y"), crs = crs)

    #--- write outputs if desired ---#
    write_samples(samples = samples, filename = filename, overwrite = overwrite)

    #--- return sf object ---#
    return(samples)
  }
}
