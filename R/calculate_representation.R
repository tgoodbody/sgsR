#' Compare sample representation within sraster strata
#'
#' @details Calculate how sraster strata are represented in existing samples
#'
#' @family calculate functions
#'
#' @inheritParams extract_strata
#' @inheritParams sample_strat
#' @param plot Logical. Plot frequency of strata coverage and sampling coverage
#' for \code{sraster} and \code{existing}. Will return a list if \code{TRUE}.
#' @param drop Numeric. Numeric value between 0-1 representing the \code{sraster} frequency
#' (\code{srasterFreq}) below which strata will be dropped from comparison (e.g..
#' This parameter can be useful for when comparing stratum where percent coverage of strata
#' may be ~ 0 percent and should be dropped. This could occur when mapping multiple stratifications.
#'
#' @return Returns a tibble where:
#' \itemize{
#' \item{strata} - \code{sraster} strata ID.
#' \item{srasterFreq} - \code{sraster} coverage frequency.
#' \item{sampleFreq} - Sampling frequency within \code{sraster} strata.
#' \item{diffFreq} - Difference between \code{srasterFreq} & \code{sampleFreq}. Positive values indicate over representation.
#' \item{nSamp} - Number of samples within each strata in \code{existing}.
#' \item{need} - \code{srasterFreq * sum(nSamp)}. Total theoretical number of required samples to be representative of strata coverage.
#' This values is rounded. It is important for the user to consider \code{diffFreq}. A small difference - e.g. 1 % -
#' in \code{sampleFreq} vs. \code{srasterFreq} frequency could make the algorithm allocate or remove samples that could likely be ignored.
#' }
#'
#' @examples
#' ### --- generate example stratification ---###
#'
#' #--- load ALS metrics from sgsR internal data ---#
#' r <- system.file("extdata", "mraster.tif", package = "sgsR")
#'
#' #--- read ALS metrics using the terra package ---#
#' mraster <- terra::rast(r)
#'
#' #--- perform stratification ---#
#' sraster <- strat_kmeans(
#'   mraster = mraster$zq90,
#'   nStrata = 6
#' )
#'
#' ### --- create existing sample network ---###
#'
#' #--- simple random sampling ---#
#' existing <- sample_srs(
#'   raster = mraster$zq90,
#'   nSamp = 100
#' )
#'
#' #--- calculate representation ---#
#'
#' calculate_representation(
#'   sraster = sraster,
#'   existing = existing
#' )
#' @author Tristan R.H. Goodbody, Martin Queinnec
#'
#' @export

calculate_representation <- function(sraster,
                                     existing,
                                     drop = NULL,
                                     plot = FALSE) {
  #--- set global vars ---#

  x <- y <- strata <- srasterFreq <- cnt <- desc <- nSamp <- sampleFreq <- diffFreq <- value <- name <- NULL

  #--- Error management ---#
  if (!inherits(sraster, "SpatRaster")) {
    stop("'sraster' must be type SpatRaster.", call. = FALSE)
  }

  suppressWarnings(
    if (!grepl("strata", names(sraster))) {
      stop("A layer name containing 'strata' does not exist within 'sraster'. Use extract_strata().", call. = FALSE)
    }
  )

  if (!inherits(existing, "data.frame") && !inherits(existing, "sf")) {
    stop("'existing' must be a data.frame or sf object.", call. = FALSE)
  }


  ### --- evaluate sample ---###
  #--- determine crs of input sraster ---#

  crs <- terra::crs(sraster, proj = TRUE)

  #--- extract covariates data from mraster ---#

  vals <- terra::as.data.frame(sraster, xy = TRUE, row.names = FALSE) %>%
    dplyr::rename(
      X = x,
      Y = y
    )

  #--- Remove NA / NaN / Inf values - calculate frequency of strata coverage ---#

  vals_mat <- vals %>%
    stats::na.omit() %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise(cnt = dplyr::n()) %>%
    dplyr::mutate(srasterFreq = round(cnt / sum(cnt), 2)) %>%
    dplyr::arrange(desc(srasterFreq))

  #--- existing ---#

  #--- avoid double strata column if existing is coming from a stratified sample function ---#
  if ("strata" %in% names(existing)) {
    existing <- existing %>%
      dplyr::select(-strata)
  }

  existing_mat <- extract_strata(sraster = sraster, existing = existing, data.frame = TRUE) %>%
    dplyr::select(strata) %>%
    stats::na.omit() %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise(nSamp = dplyr::n()) %>%
    dplyr::mutate(sampleFreq = round(nSamp / sum(nSamp), 2)) %>%
    dplyr::arrange(desc(sampleFreq))

  #--- compare ---#
  rep <- dplyr::left_join(vals_mat, existing_mat, by = "strata") %>%
    replace(is.na(.), 0) %>% # if no samples are within a strata replace NA with 0
    dplyr::mutate(diffFreq = sampleFreq - srasterFreq) %>%
    dplyr::select(strata, srasterFreq, sampleFreq, diffFreq, nSamp) %>%
    dplyr::mutate(need = ceiling(srasterFreq * sum(nSamp)) - nSamp) %>%
    dplyr::arrange(strata)

  #--- drop srasterFreq values below a certain frequency ---#
  if (!is.null(drop)) {
    if (!is.numeric(drop)) {
      stop("'drop' must be type numeric.", call. = FALSE)
    }

    if (drop > 1 | drop < 0) {
      stop("'drop' must be a numeric value between 0-1.", call. = FALSE)
    }

    rep <- rep %>%
      dplyr::filter(srasterFreq > drop)
  }

  #--- present barchart if desired ---#
  if (isTRUE(plot)) {
    p <- rep %>%
      dplyr::select(strata, sraster = srasterFreq, samples = sampleFreq) %>%
      tidyr::pivot_longer(c(2, 3)) %>%
      ggplot2::ggplot(ggplot2::aes(x = as.factor(strata), y = value, fill = name)) +
      ggplot2::geom_bar(position = "dodge", stat = "identity") +
      ggplot2::scale_fill_manual(values = c("#141414", "#5c5c5c")) +
      ggplot2::labs(
        x = "Strata",
        y = "Frequency",
        title = "Sample representation by strata",
        subtitle = "Strata coverage frequency vs. sampling frequency within strata"
      ) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.title = ggplot2::element_blank()
      )

    print(p)
  }
  return(rep)
}
