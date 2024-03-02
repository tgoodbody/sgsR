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
                                overwrite = FALSE) {
  samples <- existing %>%
    dplyr::slice_sample(., n = nSamp)

  #--- write outputs if desired ---#
  write_samples(samples = samples, filename = filename, overwrite = overwrite)

  return(samples)
}
