#' Sample allocation type and count
#'
#' @description Determine how many samples to allocate within strata.
#'
#' @family calculate functions
#'
#' @inheritParams calculate_allocation
#' @inheritParams extract_strata
#'
#' @return Returns a data.frame of:
#' \itemize{
#' \item{strata} - Strata ID.
#' \item{total} - Number of samples to be allocated. Values correspond to under representation (samples needed; positive value) or over representation
#' (too many samples; negative value) based on the \code{nSamp} provided.
#' \item{need} - Required samples per strata based on allocation method. Rounded.
#' }
#'
#' @references
#' Gregoire, T.G., & Valentine, H.T. (2007). Sampling Strategies for Natural Resources and the Environment (1st ed.).
#'  Chapman and Hall/CRC. https://doi.org/10.1201/9780203498880
#'
#' @author Tristan R.H. Goodbody
#'
#' @keywords internal
## TODO

calculate_allocation_existing <- function(existing,
                                          nSamp,
                                          allocation = "prop",
                                          weights = NULL,
                                          metric = NULL,
                                          force = FALSE) {
  if (any(!c("strata") %in% names(existing))) {
    stop("'existing must have a layer named 'strata'.", call. = FALSE)
  }

  if (!is.numeric(nSamp)) {
    stop("'nSamp' must be type numeric", call. = FALSE)
  }

  if (!any(allocation == c("prop", "optim", "equal", "manual"))) {
    stop(paste0("Unknown allocation type: '", allocation, "' provided. Please use 'prop' (default), 'optim', 'equal', or 'manual'."), call. = FALSE)
  }

  if (!is.logical(force)) {
    stop("'force' must be type logical.", call. = FALSE)
  }

  #--- set global vars ---#

  strata <- total <- eTotal <- NULL

  #--- determine which allocation algorithm to use ---#

  if (allocation != "equal") {
    #--- proportional allocation ---#

    if (allocation == "prop") {
      #--- error handling when allocation algorithm is 'prop' ---#

      if (!is.null(metric)) {
        message("'metric' was specified but 'allocation = prop' - did you mean to use 'allocation = optim'?")
      }

      if (!is.null(weights)) {
        message("'weights' was specified but 'allocation = prop' - did you mean to use 'allocation = manual'?")
      }

      toSample <- allocate_existing_prop(existing = existing, nSamp = nSamp)
    }

    #--- optimal allocation ---#

    if (allocation == "optim") {
      #--- error handling when allocation algorithm is 'optim' ---#

      if (!is.null(weights)) {
        message("'weights' was specified but 'allocation = optim' - did you mean to use 'allocation = manual'?")
      }

      toSample <- allocate_existing_optim(existing = existing, metric = metric, nSamp = nSamp)
    }

    #--- manual allocation ---#
    if (allocation == "manual") {
      #--- error handling when allocation algorithm is 'manual' ---#

      if (!is.null(metric)) {
        message("'metric' was specified but 'allocation = manual' - did you mean to use 'allocation = optim'?")
      }

      toSample <- allocate_existing_manual(existing = existing, nSamp = nSamp, weights = weights)
    }

    #--- calculate total samples allocated ---#

    tot <- sum(toSample$total)

    if (isFALSE(force) && tot != nSamp) {
      message(paste0("nSamp of ", nSamp, " is not perfectly divisible based on strata distribution. nSamp of ", tot, " will be returned. Use 'force = TRUE' to brute force to ", nSamp, "."))
    }
  } else {
    #--- equal allocation ---#

    if (!is.null(weights)) {
      message("'weights' was specified but 'allocation = equal' - did you mean to use 'allocation = manual'?")
    }

    if (!is.null(metric)) {
      message("'metric' was specified but 'allocation = equal' - did you mean to use 'allocation = optim'?")
    }

    toSample <- allocate_existing_equal(existing = existing, nSamp = nSamp)

    tot <- unique(toSample$total)
  }

  #--- determine whether there is a difference between 'nSamp' and the number of allocated samples with each stratum ---#

  diff <- tot - nSamp

  if (allocation != "equal") {
    if (diff != 0) {
      #--- adjust sample count to force the user defined number ---#

      if (isTRUE(force)) {
        message(paste0("Forcing ", nSamp, " total samples."))

        #--- if samples need to be removed ---#

        toSample <- allocate_force(toSample = toSample, nSamp = nSamp, diff = diff)
      }
    }
  } else {
    if (force == TRUE) {
      message("`force = TRUE` has no effect when `allocation = equal'. Ignoring.")
    }
  }

  toSample
}
