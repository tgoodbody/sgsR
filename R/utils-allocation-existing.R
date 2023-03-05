#' Allocating strata: `existing`
#'
#' Allocation algorithms based on `existing`
#'
#' @inheritParams calculate_allocation_existing
#' @family allocation
#' @name allocating
#' @return Data frame of allocated samples by stratum. Used internally within \code{sample_existing(type = "strat")}.
NULL

#' @export
#' @rdname allocating
#' @family allocation
#' @keywords internal

allocate_existing_prop <- function(existing,
                                   nSamp) {
  #--- define global vars ---#

  strata <- count <- freq <- total <- NULL

  message("Implementing proportional allocation of samples.")

  #--- generate vals data.frame ---#

  existing %>%
    sf::st_drop_geometry() %>%
    dplyr::select(strata) %>%
    stats::na.omit() %>%
    dplyr::group_by(strata) %>%
    dplyr::summarize(count = dplyr::n()) %>%
    dplyr::mutate(
      freq = count / sum(count),
      total = freq * nSamp
    ) %>%
    dplyr::mutate(total = replace(total, total < 1, 1)) %>%
    dplyr::mutate(total = round(total, digits = 0)) %>%
    dplyr::select(strata, total)
}

#' @export
#' @rdname allocating
#' @family allocation
#' @keywords internal

allocate_existing_optim <- function(existing,
                                    metric,
                                    nSamp) {
  #--- define global vars ---#

  v_sd <- strata <- count <- total <- denom <- NULL

  #--- error handling when allocation algorithm is 'optim' ---#

  if (is.null(metric)) {
    stop("'metric' must be supplied if 'allocation = optim'.", call. = FALSE)
  }

  #--- if there is only 1 band in mraster use it as default ---#

  if (length(metric) > 1) {
    stop("Multiple character strings detected in 'metric'. Please define a singular metric for allocation.", call. = FALSE)
  }

  if (!metric %in% names(existing)) {
    stop(paste0("No column named ", metric, " in 'existing'.", call. = FALSE))
  }

  message(paste0("Implementing optimal allocation of samples based on variability of '", metric, "'."))

  #--- merge sraster and mraster together ---#

  existing %>%
    sf::st_drop_geometry() %>%
    dplyr::select(strata, {{ metric }}) %>%
    stats::na.omit() %>%
    dplyr::group_by(strata) %>%
    dplyr::summarize(
      v_sd = sd(.[[{{ metric }}]]),
      count = dplyr::n()
    ) %>%
    dplyr::mutate(denom = sum(count * v_sd)) %>%
    dplyr::rowwise() %>%
    #--- optimal allocation (equal sampling cost) equation. See Gregoire & Valentine (2007) Section 5.4.4 ---#
    dplyr::mutate(total = round(nSamp * ((count * v_sd) / denom), digits = 0)) %>%
    dplyr::select(strata, total)
}

#' @export
#' @rdname allocating
#' @family allocation
#' @keywords internal

allocate_existing_manual <- function(existing,
                                     nSamp,
                                     weights) {
  #--- define global vars ---#

  strata <- total <- NULL

  #--- error handling when allocation algorithm is 'manual' ---#

  if (is.null(weights)) {
    stop("'weights' must be defined if 'allocation = manual'.", call. = FALSE)
  }

  if (!is.numeric(weights)) {
    stop("'weights' must be a numeric vector.", call. = FALSE)
  }

  if (sum(weights) != 1) {
    stop("'weights' must add up to 1.", call. = FALSE)
  }

  message("Implementing allocation of samples based on user-defined weights.")

  #--- generate vals data.frame ---#

  vals <- existing %>%
    sf::st_drop_geometry() %>%
    stats::na.omit()

  if (length(weights) != length(unique(vals$strata))) {
    stop("'weights' must be the same length as the number of strata in 'sraster'.", call. = FALSE)
  }

  #--- determine number of samples within each strata ---#
  vals %>%
    dplyr::group_by(strata) %>%
    dplyr::summarize(count = dplyr::n()) %>%
    dplyr::mutate(
      weights = weights,
      total = nSamp * weights
    ) %>%
    dplyr::mutate(total = replace(total, total < 1, 1)) %>%
    dplyr::mutate(total = round(total, digits = 0)) %>%
    dplyr::select(strata, total)
}

#' @export
#' @rdname allocating
#' @family allocation
#' @keywords internal

allocate_existing_equal <- function(existing,
                                    nSamp) {
  message("Implementing equal allocation of samples.")

  #--- define global vars ---#

  strata <- NULL

  #--- generate vals data.frame ---#

  toSample <- existing %>%
    sf::st_drop_geometry() %>%
    stats::na.omit() %>%
    dplyr::group_by(strata) %>%
    dplyr::summarize(total = nSamp)

  toSample
}
