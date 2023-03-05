#' Randomly sample from an existing dataset
#'
#' This function takes an existing \code{sf} object and returns a random sub-sample of size \code{nSamp}.
#'
#' @inheritParams extract_strata
#' @inheritParams sample_srs
#' 
#' @return A data frame or spatial data frame containing the random sample.
#'
#' @keywords internal
sample_existing_srs <- function(existing,
                                nSamp,
                                filename = NULL,
                                overwrite = NULL) {
  samples <- existing %>%
    dplyr::slice_sample(., n = nSamp)

  #--- write to disc ---#

  if (!is.null(filename)) {
    if (!is.character(filename)) {
      stop("'filename' must be a file path character string.", call. = FALSE)
    }

    if (!is.logical(overwrite)) {
      stop("'overwrite' must be type logical.", call. = FALSE)
    }

    if (file.exists(filename) & isFALSE(overwrite)) {
      stop(paste0("'", filename, "' already exists and overwrite = FALSE."), call. = FALSE)
    }

    sf::st_write(samples, filename, delete_layer = overwrite)
    message("Output samples written to disc.")
  }

  return(samples)
}
