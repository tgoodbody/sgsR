#' Quantiles stratification
#'
#' @description Stratify metric raster using metric quantiles.
#' @family stratify functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams strat_breaks
#' @param metric Numeric/Character. Index or name of primary covariate within \code{mraster} to stratify.
#' @param metric2 Numeric/Character. Index or name of secondary covariate within \code{mraster} to stratify.
#' @param nQuant Numeric. Number of quantiles to stratify primary covariate.
#' @param nQuant2 Numeric. Number of quantiles to stratify secondary covariate.
#' @param samp Numeric. Determines proportion of cells to plot in scatterplot - see \code{values}
#' for strata visualization. Lower values reduce visualization time.
#'
#'
#' @return Returns an output stratification \code{spatRaster} or a list when \code{details = TRUE}.
#'
#' When a list is returned:
#' \enumerate{
#' \item \code{details} is a list output of the \code{\link[stats]{prcomp}} function
#' \item \code{raster} is a stratified \code{spatRaster} based on quantiles
#' \item \code{plot} is a \code{ggplot} histogram / scatter plot object (depends on whether metric2 was supplied).
#' Histogram shows distribution and break points while scatter plot shows colour coded and strata boundaries.
#' }
#'
#' @examples
#' #--- Load raster and existing plots---#
#' r <- system.file("extdata", "wall_metrics_small.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' strat_quantiles(mraster = mr, 
#'                 metric = 4, 
#'                 nQuant = 10,
#'                 plot = TRUE, 
#'                 details = TRUE)
#' 
#' strat_quantiles(mraster = mr, 
#'                 metric = "zsd",
#'                 metric2 = "zq95", 
#'                 nQuant = 3, 
#'                 nQuant2 = 4)
#' 
#' strat_quantiles(mraster = mr, 
#'                 metric = 1, 
#'                 metric2 = "zsd",
#'                 nQuant = 2, 
#'                 nQuant2 = 2, 
#'                 filename = tempfile(fileext = ".tif"))
#'                 
#' @author Tristan R.H. Goodbody
#'
#' @export

strat_quantiles <- function(mraster,
                            metric = NULL,
                            metric2 = NULL,
                            nQuant,
                            nQuant2 = NULL,
                            plot = FALSE,
                            details = FALSE,
                            samp = 1,
                            filename = NULL,
                            overwrite = FALSE,
                            ...) {

  #--- Set global vars ---#

  class1 <- class2 <- NULL

  #--- error handling ---#
  if (!inherits(mraster, "SpatRaster")) {
    stop("all specified bands must be type SpatRaster", call. = FALSE)
  }

  if (!is.numeric(nQuant)) {
    stop("'nQuant' must be type numeric")
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical")
  }

  if (!is.numeric(samp)) {
    stop("'samp' must be type numeric")
  }

  if (!is.logical(details)) {
    stop("'details' must be type logical")
  }

  if (is.null(metric2)) {
    if (!is.null(nQuant2)) {
      message("You are stratifying with only 1 metric but specified 'nQuant2' - ignoring.")
    }

    #--- if there is only 1 metric in the raster use it as default ---#

    if (terra::nlyr(mraster) == 1) {

      #--- Extract values from mraster ---#

      vals <- terra::values(mraster)

      #--- set name of raster band to 'metric' ---#

      metric <- names(mraster)
    } else {

      #--- subset metric based on whether it is a character of number ---#

      if (is.null(metric)) {
        stop(" multiple layers detected in 'mraster'. Please define a 'metric' to stratify")
      } else {

        #--- Numeric ---#

        if (is.numeric(metric)) {
          if ((metric) > (terra::nlyr(mraster)) | metric < 0) {
            stop("'metric' index doest not exist within 'mraster'")
          }

          #--- Character ---#
        } else if (is.character(metric)) {
          if (!metric %in% names(mraster)) {
            stop(glue::glue("'mraster' must have an attribute named {metric}."))
          }

          metric <- which(names(vals) == metric)
        }
      }

      #--- Extract values from mraster ---#

      vals <- terra::subset(mraster, metric) %>%
        terra::values()
    }

    vals[!is.finite(vals)] <- NA

    #--- Determine index of each cell so to map values correctly without NA ---#

    idx <- !is.na(vals)

    #--- Remove NA / NaN / Inf values ---#

    df <- vals %>%
      as.data.frame() %>%
      dplyr::filter(!is.na(.))

    metric <- as.character(names(df))

    metric <- ggplot2::ensym(metric)

    #--- Split metric distribution in to number specified by 'breaks' ---#

    dfc <- df %>%
      dplyr::mutate(class = dplyr::ntile(!!metric, nQuant))

    #--- convert back to original mraster extent ---#

    vals[idx] <- dfc$class

    #--- set newly stratified values ---#

    rout <- terra::setValues(mraster[[1]], vals)
    names(rout) <- "strata"
  } else {

    #--- subset metric2 based on whether it is a character of number ---#

    if (is.null(nQuant2)) {
      stop("If using 2 metrics to stratify, 'nQuant2' must be defined")
    }

    #--- metric 1 ---#

    #--- numeric ---#
    if (is.numeric(metric)) {
      if ((metric) > (terra::nlyr(mraster)) | metric < 0) {
        stop("'metric' index doest not exist within 'mraster'")
      }

      #--- Character ---#
    } else if (is.character(metric)) {
      if (!metric %in% names(mraster)) {
        stop(glue::glue("'mraster' must have an attribute named {metric}."))
      }

      metric <- which(names(mraster) == metric)
    }

    #--- metric 2 ---#

    #--- Numeric ---#

    if (is.numeric(metric2)) {
      if ((metric2) > (terra::nlyr(mraster)) | metric2 < 0) {
        stop("'metric' index doest not exist within 'mraster'")
      }

      #--- Character ---#
    } else if (is.character(metric2)) {
      if (!metric2 %in% names(mraster)) {
        stop(glue::glue("'mraster' must have an attribute named {metric}."))
      }

      metric2 <- which(names(mraster) == metric2)
    }
    #--- Extract values from mraster ---#

    vals <- terra::subset(mraster, c(metric, metric2)) %>%
      terra::values()

    vals[!is.finite(vals)] <- NA

    #--- Determine index of each cell so to map values correctly without NA ---#

    idx <- is.finite(vals[, 1]) & is.finite(vals[, 2])

    #--- Remove NA / NaN / Inf values ---#

    df <- vals %>%
      as.data.frame() %>%
      dplyr::filter(!is.na(.))

    nm <- as.character(names(df))

    metric <- nm[1]
    metric2 <- nm[2]

    metric <- ggplot2::ensym(metric)
    metric2 <- ggplot2::ensym(metric2)

    #--- Split metric distribution in to number specified by 'breaks' ---#

    dfc <- df %>%
      #--- define nQuant classes ---#
      dplyr::mutate(class1 = dplyr::ntile(!!metric, nQuant)) %>%
      #--- group by class to sub stratify ---#
      dplyr::group_by(class1) %>%
      #--- define nQuant2 classes ---#
      dplyr::mutate(class2 = dplyr::ntile(!!metric2, nQuant2)) %>%
      #--- combine classes ---#
      dplyr::group_by(class1, class2) %>%
      #--- establish newly formed unique class ---#
      dplyr::mutate(class = dplyr::cur_group_id())

    #--- convert back to original mraster extent ---#

    suppressWarnings(vals[, 1][idx] <- dfc$class)

    #--- set newly stratified values ---#

    rout <- terra::setValues(mraster[[1]], vals[, 1])
    names(rout) <- "strata"
  }

  if (isTRUE(plot)) {
    if (samp > 1 | samp < 0) {
      stop("'samp' must be between 0 and 1")
    }

    if (is.null(metric2)) {

      #--- output histogram with quantile breaks ---#

      #--- determine numeric breaks ---#

      breaks <- dfc %>%
        dplyr::group_by(class) %>%
        dplyr::summarize(breaks = max(!!metric)) %>%
        dplyr::select(breaks) %>%
        as.data.frame()

      breaks <- breaks[1:(nQuant - 1), ]

      df.p <- dfc %>%
        dplyr::select(metric)

      #--- plot histogram of metric with associated break lines ---#

      p <- ggplot2::ggplot(df.p, ggplot2::aes(!!metric)) +
        ggplot2::geom_histogram() +
        ggplot2::geom_vline(xintercept = breaks, linetype = "dashed") +
        ggplot2::ggtitle(glue::glue('{metric} histogram with defined breaks'))

      print(p)

      terra::plot(rout, main = "Classes")
    } else {
      terra::plot(rout, main = "Classes")
    }

    #--- write file to disc ---#

    if (!is.null(filename)) {
      terra::writeRaster(rout, filename, overwrite = overwrite, ...)
    }

    #--- Output based on 'details' to return raster alone or list with details ---#

    if (isTRUE(details)) {
      if (!is.null(metric2)) {

        #--- create classplot summary ---#

        coordsgrps <- dfc %>%
          dplyr::group_by(class) %>%
          dplyr::arrange(class) %>%
          tidyr::nest() %>%
          dplyr::ungroup()

        p <- classPlot(
          dfc = dfc,
          coordsgrps = coordsgrps,
          metric = metric,
          metric2 = metric2,
          samp = samp
        )
      }

      #--- output metrics details along with stratification raster ---#

      out <- list(
        details = dfc,
        raster = rout,
        plot = p
      )

      return(out)
    } else {

      #--- just output raster ---#

      return(rout)
    }
  }

  return(rout)
}
