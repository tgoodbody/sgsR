#' Sample Existing Data Based on Strata
#'
#' This function takes a data frame of existing data, a data frame of desired sample sizes
#' for each strata, the number of samples to take, and optionally a file name and overwrite
#' parameter. It returns a sample of the existing data based on the desired sample sizes for
#' each strata, and optionally writes the resulting samples to a file.
#'
#' @inheritParams extract_strata
#' @inheritParams sample_strat
#' @param toSample A data frame specifying the desired sample sizes for each strata.
#'
#' @return An \code{sf} object that is a sub-sample of \code{existing}
#'
#' @keywords internal
#' @noRd
sample_existing_strat <- function(existing,
                                  toSample,
                                  nSamp,
                                  filename = NULL,
                                  overwrite = FALSE) {
  strata <- unique(existing$strata)


  samples <- mapply(take_samples, strata = strata, MoreArgs = list(existing, toSample), SIMPLIFY = FALSE) %>%
    dplyr::bind_rows()


  #--- write outputs if desired ---#
  write_samples(samples = samples, filename = filename, overwrite = overwrite)

  return(samples)
}

#' Take Samples Based on Strata
#'
#' This function takes a data frame of existing data, a data frame of desired sample sizes
#' for each strata, and a strata variable name, and returns a sample of the existing data
#' based on the sample sizes for the specified strata.
#'
#' @param existing A data frame containing existing data.
#' @param toSample A data frame specifying the desired sample sizes for each strata.
#' @param strata A string specifying the name of the variable used to define strata.
#'
#' @return A data frame containing a sample of the existing data based on the sample sizes
#' for the specified strata.
#'
#' @keywords internal
#' @noRd
take_samples <- function(existing, toSample, strata) {
  total <- NULL

  toTake <- toSample %>%
    dplyr::filter(strata == {{ strata }}) %>%
    dplyr::select(total) %>%
    dplyr::pull()

  tryCatch(
    {
      existing %>%
        dplyr::filter(strata == {{ strata }}) %>%
        dplyr::slice_sample(., n = toTake)
    },
    error = function(e) {
      if (grepl("cannot take a sample larger than the population", e$message)) {
        stop("Error: ", e$message)
      } else {
        stop(paste0("Not enough samples in strata: ", strata, " to take: ", toTake, " sample units."), call. = FALSE)
      }
    }
  )
}
