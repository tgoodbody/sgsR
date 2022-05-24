#' Sample size determination
#'
#' @description Determine a samples size for simple random sampling using relative standard error
#'
#' @family calculate functions
#'
#' @param mraster spatRaster. Metrics raster. All values must be numeric.
#' @param rse Numeric. Desired relative standard error (coefficient of variation of the mean)
#' threshold to determine sample size.
#' @param start Numeric. First rse value to begin rse sequence. default = \code{0.01}.
#' @param end Numeric. Final rse value to end rse sequence. default = \code{0.05}.
#' @param increment Numeric. Value to increment between \code{start} & \code{end}. default = \code{0.001}.
#' @param plot Logical. if \code{TRUE} output graphical representation of estimated sample size vs. rse.
#'
#' @return A data.frame of sample size and rse by raster variable.
#'
#' @examples
#'
#' #--- Load raster ---#
#' r <- system.file("extdata", "wall_metrics.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' calculate_sampsize(
#'   mraster = mr,
#'   rse = 0.01
#' )
#'
#' calculate_sampsize(
#'   mraster = mr
#' )
#'
#' calculate_sampsize(
#'   mraster = mr,
#'   rse = 0.025,
#'   start = 0.01,
#'   end = 0.08,
#'   increment = 0.01
#' )
#'
#' #--- higher variance leads to need for more samples ---#
#' @author Tristan R.H. Goodbody
#'
#' @note
#'
#' \deqn{rse = (100 * SE) / mean)}
#'
#' Where:
#'
#' \itemize{
#'  \item{SE}{ - Standard error of the mean}
#'  \item{s}{ - Standard deviation of the observations}
#'  \item{n}{ - Number of observations}
#' }
#'
#' @references
#' Benedetti, R., Piersimoni, F., & Postiglione, P. (2015).
#' Sampling spatial units for agricultural surveys. pp 202-203. Berlin: Springer.
#'
#' @export


calculate_sampsize <- function(mraster,
                               rse = NULL,
                               start = 0.01,
                               end = 0.05,
                               increment = 0.001,
                               plot = FALSE) {

  #--- set global vars ---#

  rse_var <- nSamp <- var <- rse_var_dif <- NULL

  #--- error management ---#

  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster")
  }

  if (!is.numeric(start)) {
    stop("'start' must be type numeric")
  }

  if (!is.numeric(end)) {
    stop("'end' must be type numeric")
  }

  if (!is.numeric(increment)) {
    stop("'increment' must be type increment")
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical")
  }

  #--- convert raster to data.frame ---#
  vals <- terra::values(mraster, dataframe = TRUE) %>%
    na.omit()

  if (!any(apply(vals, 2, is.numeric))) {
    stop("'mraster' must contain all numeric values.")
  }

  #--- calculate adequate sample sizes based on sequence of relative standard errors ---#

  rse_seq <- seq(start, end, increment)

  sampsize <- apply(X = vals, MARGIN = 2, FUN = size_calculation, N = nrow(vals), rse = rse_seq)

  #--- apply appropriate variable indices ---#

  nm <- names(sampsize)

  for (n in 1:length(nm)) {
    sampsize[[n]]$var <- nm[n]
  }

  #--- clean data ---#

  sampsize <- do.call(rbind, sampsize)
  rownames(sampsize) <- c()

  if (is.null(rse)) {

    #--- plot ---#

    if (isTRUE(plot)) {
      p <- ggplot2::ggplot(sampsize, ggplot2::aes(x = rse_var, y = nSamp)) +
        ggplot2::geom_point(colour = "#333333") +
        ggplot2::geom_line() +
        ggplot2::facet_wrap(. ~ var, scales = "free") +
        ggplot2::ggtitle("Samples sizes by rse") +
        ggplot2::xlab("Required relative standard error") +
        ggplot2::ylab("Sample size") +
        ggplot2::theme_bw()
    }

    names(sampsize) <- c("nSamp", "rse", "var")

    #--- change name for single output ---#

    lines <- sampsize
  } else {

    #--- checks for rse ---#

    if (!is.numeric(rse)) {
      stop("'rse' must be type numeric")
    }

    if (rse < 0) {
      stop("`rse` must be > 0")
    }

    if (rse > 0.15) {
      message("`rse` > 0.15 -  are you sure you want this?")
    }

    if (rse < start | rse > end) {
      stop("'rse' must be >= `start` and <= `end`")
    }

    #--- determine if increment and rse are divisble ---#

    if (rse %% increment != 0) {
      lines <- sampsize %>%
        dplyr::group_by(var) %>%
        dplyr::mutate(rse_var_dif = abs(rse - rse_var)) %>%
        dplyr::filter(rse_var_dif == min(rse_var_dif)) %>%
        dplyr::select(-rse_var_dif)

      message(paste0("'rse' not perfectly divisible by 'incremenent. Selecting closest sample size (rse = ", unique(lines$rse_var), ") based on values."))
    } else {
      lines <- sampsize %>%
        dplyr::filter(rse_var == rse)
    }

    #--- plot ---#

    if (isTRUE(plot)) {
      p <- ggplot2::ggplot(sampsize, ggplot2::aes(x = rse_var, y = nSamp)) +
        ggplot2::geom_point(colour = "#333333") +
        ggplot2::geom_line() +
        ggplot2::geom_segment(data = lines, ggplot2::aes(x = -Inf, xend = rse_var, y = nSamp, yend = nSamp), linetype = "dashed", colour = "red") +
        ggplot2::geom_segment(data = lines, ggplot2::aes(x = rse_var, xend = rse_var, y = -Inf, yend = nSamp), linetype = "dashed", colour = "red") +
        ggplot2::geom_text(data = lines, mapping = ggplot2::aes(label = paste0("nSamp = ", nSamp), x = Inf, y = Inf, vjust = 2, hjust = 1.2), colour = "red") +
        ggplot2::xlim(min(sampsize$rse_var), max(sampsize$rse_var)) +
        ggplot2::facet_wrap(. ~ var, scales = "free") +
        ggplot2::ggtitle(paste0("Samples size with rse = ", unique(lines$rse_var))) +
        ggplot2::xlab("Required relative standard error") +
        ggplot2::ylab("Sample size") +
        ggplot2::theme_bw()
    }

    names(lines) <- c("nSamp", "rse", "var")
  }

  if (exists("p")) {
    out <- list(
      nSamp = lines,
      plot = p
    )
  } else {
    out <- lines
  }

  return(out)
}

#--- determine sample size based on rse value ---#

size_calculation <- function(N,
                             mvals,
                             rse) {

  #--- determine adequate samples size based on the relative standard error defined ---#

  nSamp <- ceiling((N^2 * var(mvals)) / (rse^2 * sum(mvals)^2 + N * var(mvals)))

  out <- data.frame(nSamp = nSamp, rse_var = rse)

  out
}
