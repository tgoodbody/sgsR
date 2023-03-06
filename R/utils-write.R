#' Write
#'
#' @inheritParams sample_systematic
#' @param raster SpatRaster. Raster to be written to disc
#' @family write
#' @name write
#' @return No explicit output but files are written to disc if specified.
NULL

#' Write samples
#' @family write
#' @rdname write
#' @keywords internal

write_samples <- function(samples, filename = NULL, overwrite = FALSE) {
  if (!is.null(filename)) {
    if (!is.character(filename)) {
      stop("'filename' must be a file path character string.", call. = FALSE)
    }

    if (!is.logical(overwrite)) {
      stop("'overwrite' must be type logical.", call. = FALSE)
    }

    if (file.exists(filename) & isFALSE(overwrite)) {
      stop(paste0("'filename' already exists and overwrite = FALSE."), call. = FALSE)
    }

    sf::st_write(samples, filename, delete_layer = overwrite)
    message("Output samples written to disc.")
  }
}

#' Write table
#' @family write
#' @rdname write
#' @keywords internal

write_samples_df <- function(samples, filename = NULL, overwrite = FALSE) {
  if (!is.null(filename)) {
    if (!is.character(filename)) {
      stop("'filename' must be type character.", call. = FALSE)
    }

    if (!is.logical(overwrite)) {
      stop("'overwrite' must be type logical.", call. = FALSE)
    }

    #--- append and overwrite are opposites .. need to invert them for csv writing ---#

    if (file.exists(filename) & isFALSE(overwrite)) {
      stop(paste0("'filename' already exists and overwrite = FALSE."), call. = FALSE)
    }

    utils::write.table(x = samples, file = filename, append = !overwrite)
    message("Output samples written to disc.")
  }
}

#' Write raster
#' @family write
#' @rdname write
#' @keywords internal

write_raster <- function(raster, filename = NULL, overwrite = FALSE) {
  if (!is.null(filename)) {
    if (!is.character(filename)) {
      stop("'filename' must be a file path character string.", call. = FALSE)
    }

    if (!is.logical(overwrite)) {
      stop("'overwrite' must be type logical.", call. = FALSE)
    }

    if (file.exists(filename) & isFALSE(overwrite)) {
      stop(paste0("'filename' already exists and overwrite = FALSE."), call. = FALSE)
    }

    terra::writeRaster(x = raster, filename = filename, overwrite = overwrite)
    message("Output raster written to disc.")
  }
}
