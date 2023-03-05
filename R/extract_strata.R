#' Extract strata
#'
#' @description Extract stratum values to existing samples
#'
#' @family extract functions
#'
#' @inheritParams sample_systematic
#'
#' @param sraster spatRaster. Stratification raster.
#' @param existing sf 'POINT'.  Existing plot network.
#' @param data.frame Logical. Output as data.frame if \code{TRUE}
#' @param quiet Logical. If \code{TRUE} the user will not get messages
#' about samples with \code{NA} values.
#'
#' @return An sf or data.frame object of samples with strata attribute.
#'
#' @note
#'
#' If \code{data.frame = TRUE} output will be written using \code{\link[utils]{write.table}}
#'
#' @examples
#' #--- Load sraster ---#
#' r <- system.file("extdata", "sraster.tif", package = "sgsR")
#' sr <- terra::rast(r)
#'
#' #--- load existing samples ---#
#' e <- system.file("extdata", "existing.shp", package = "sgsR")
#' e <- sf::st_read(e)
#'
#' extract_strata(
#'   sraster = sr,
#'   existing = e
#' )
#'
#' @author Tristan R.H. Goodbody
#'
#' @export

extract_strata <- function(sraster,
                           existing,
                           quiet = FALSE,
                           data.frame = FALSE,
                           filename = NULL,
                           overwrite = FALSE) {
  #--- Set global vars ---#

  x <- y <- X <- Y <- strata <- ID <- geometry <- NULL

  #--- Error management ---#

  if (!inherits(sraster, "SpatRaster")) {
    stop("'sraster' must be type SpatRaster.", call. = FALSE)
  }

  if (terra::nlyr(sraster) > 1) {
    stop("'sraster' must have a single layer named 'strata'.", call. = FALSE)
  }

  if (any(!c("strata") %in% names(sraster))) {
    stop("'sraster' must have a layer named 'strata'.", call. = FALSE)
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
    crs <- terra::crs(sraster)

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

  vals <- terra::extract(sraster, xy)

  #--- when dataframe is input "ID" is appended to vals -- remove it ---#

  if ("ID" %in% names(vals)) {
    vals <- vals %>% dplyr::select(-ID)
  }

  #--- check that extractions has produced some values -- if not tell the user ---#

  if (all(!complete.cases(vals))) {
    stop("'existing' only extracts NA values. Ensure that 'existing' overlaps with 'sraster'.", call. = FALSE)
  }

  #--- if existing samples are not linked with a stratum ---#
  if (any(!complete.cases(vals))) {
    nNA <- vals %>%
      dplyr::filter(!complete.cases(strata)) %>%
      dplyr::tally() %>%
      dplyr::pull()

    if (isFALSE(quiet)) {
      message(paste0(nNA, " samples are located where strata values are NA."))
    }
  }

  #--- output either data.frame or sf object ---#

  if (isTRUE(data.frame)) {
    samples <- cbind(xy, vals, existing)

    if ("geometry" %in% names(samples)) {
      samples <- samples %>%
        dplyr::select(-geometry)
    }

    if (!is.null(filename)) {
      if (!is.character(filename)) {
        stop("'filename' must be type character.", call. = FALSE)
      }

      if (!is.logical(overwrite)) {
        stop("'overwrite' must be type logical.", call. = FALSE)
      }

      #--- append and overwrite are opposites .. need to invert them for csv writing ---#

      if (file.exists(filename) & isFALSE(overwrite)) {
        stop(paste0("'", filename, "' already exists and overwrite = FALSE"))
      }

      utils::write.table(x = samples, file = filename, append = !overwrite)
      message("Output samples written to disc.")
    }

    #--- return data.frame ---#
    return(samples)
  } else {
    #--- convert coordinates to a sf object ---#
    samples <- cbind(xy, vals, existing) %>%
      sf::st_as_sf(., coords = c("X", "Y"), crs = crs)

    if (!is.null(filename)) {
      if (!is.character(filename)) {
        stop("'filename' must be type character.", call. = FALSE)
      }

      if (!is.logical(overwrite)) {
        stop("'overwrite' must be type logical.", call. = FALSE)
      }

      if (file.exists(filename) & isFALSE(overwrite)) {
        stop(paste0("'", filename, "' already exists and overwrite = FALSE"))
      }

      sf::st_write(samples, filename, delete_layer = overwrite)
      message("Output samples written to disc.")
    }

    #--- return sf object ---#
    return(samples)
  }
}
