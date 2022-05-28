#' Breaks stratification
#'
#' @description Stratify metrics raster using user defined breaks
#' @family stratify functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams sample_systematic
#' @param mraster Spatraster. Primary covariate raster to stratify.
#' @param mraster2 Spatraster. Secondary covariate raster to stratify.
#' @param breaks Numeric. Vector of breakpoints for \code{mraster}.
#' @param breaks2 Numeric. Vector of breakpoints for \code{mraster2} (if provided).
#' @param filename Character. Path to write stratified raster to disc.
#' @param overwrite Logical. Specify whether \code{filename} should be overwritten on disc.
#' @param ... Additional arguments for writing files. See \code{\link[terra]{writeRaster}}.
#'
#' @return Returns an output stratification \code{spatRaster} or a list when \code{details = TRUE}.
#'
#' When a list is returned:
#' \enumerate{
#' \item \code{details} is a list output of the \code{\link[stats]{prcomp}} function
#' \item \code{raster} is a stratified \code{spatRaster} based on quantiles
#' \item \code{plot} is a \code{ggplot} histogram object showing distribution and break points.
#' }
#'
#' @examples
#' #--- Load raster ---#
#' r <- system.file("extdata", "mraster.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' #--- create vector breaks ---#
#' br.max <- c(3, 5, 11, 18)
#' br.sd <- c(1, 2, 5)
#'
#' strat_breaks(
#'   mraster = mr$zq90,
#'   breaks = br.max,
#'   plot = TRUE,
#'   details = TRUE
#' )
#'
#' strat_breaks(
#'   mraster = mr$zq90,
#'   mraster2 = mr$zsd,
#'   breaks = br.max,
#'   breaks2 = br.sd
#' )
#' @author Tristan R.H. Goodbody
#'
#' @export

strat_breaks <- function(mraster,
                         mraster2 = NULL,
                         breaks,
                         breaks2 = NULL,
                         plot = FALSE,
                         details = FALSE,
                         filename = NULL,
                         overwrite = FALSE,
                         ...) {

  #--- Set global vars ---#
  from <- strata <- strata2 <- val <- brk <- NULL

  #--- Error management ---#

  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster")
  }

  if (!is.numeric(breaks)) {
    stop("'breaks' must be type numeric")
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical")
  }

  if (!is.logical(details)) {
    stop("'details' must be type logical")
  }

  #--- if there is only 1 band in mraster use it as default ---#

  if (terra::nlyr(mraster) == 1) {
    rastermetric <- mraster
  } else {
    stop("Multiple layers detected in 'mraster'. Please define a singular band to stratify.")
  }

  if (!is.null(mraster2)) {
    if (!inherits(mraster2, "SpatRaster")) {
      stop("'mraster2' must be type SpatRaster")
    }

    if (!all.equal(terra::ext(mraster), terra::ext(mraster2))) {
      stop("Extents of 'mraster' and 'mraster2' do not match.")
    }

    if (!all.equal(terra::res(mraster), terra::res(mraster2))) {
      stop("Spatial resolutions of 'mraster' and 'mraster2' do not match.")
    }

    #--- subset mraster2 based on whether it is a character of number ---#

    if (is.null(breaks2)) {
      stop("If using mraster2 to stratify, 'breaks2' must be defined")
    }

    #--- if there is only 1 band in mraster2 use it as default ---#

    if (terra::nlyr(mraster2) == 1) {
      rastermetric2 <- mraster2
    } else {
      stop("Multiple layers detected in 'mraster2'. Please define a singular band to stratify.")
    }
  }

  minmax <- terra::minmax(rastermetric)

  if (any(breaks < minmax[1])) {
    message("'breaks' contains values < the minimum 'mraster' value.")
  }

  if (any(breaks > minmax[2])) {
    message("'breaks' contains values > the maximum 'mraster' value.")
  }

  #--- apply breaks to primary metric ---#

  #--- reclassify values based on breaks ---#

  breaks_m <- data.frame(from = c(-Inf, breaks, Inf)) %>%
    dplyr::mutate(
      to = dplyr::lead(from),
      becomes = seq(1:length(from))
    ) %>%
    stats::na.omit() %>%
    as.matrix()

  rcl <- terra::classify(rastermetric, breaks_m, othersNA = TRUE)
  names(rcl) <- "strata"


  #--- if secondary metric is provided ---#

  if (!is.null(mraster2)) {

    #--- Determine if breaks are > / < min and max values of metrics ---#

    minmax2 <- terra::minmax(rastermetric2)

    if (any(breaks2 < minmax2[1])) {
      stop("'breaks2' contains values < the minimum 'mraster2' value.")
    }

    if (any(breaks2 > minmax2[2])) {
      stop("'breaks2' contains values > the maximum 'mraster2' value.")
    }

    #--- reclassify values based on breaks ---#

    breaks2_m <- data.frame(from = c(-Inf, breaks2, Inf)) %>%
      dplyr::mutate(
        to = dplyr::lead(from),
        becomes = seq(1:length(from))
      ) %>%
      stats::na.omit() %>%
      as.matrix()

    rcl2 <- terra::classify(rastermetric2, breaks2_m, othersNA = TRUE)
    names(rcl2) <- "strata2"

    #--- stack rcl and rcl2

    rstack <- c(rcl, rcl2)

    breaks_c <- terra::as.data.frame(rstack, xy = TRUE)

    breaks_c <- breaks_c %>%
      dplyr::group_by(strata, strata2) %>%
      #--- establish newly formed unique class ---#
      dplyr::mutate(class = dplyr::cur_group_id()) %>%
      dplyr::ungroup()

    idx <- terra::cellFromXY(rcl, cbind(breaks_c$x, breaks_c$y))

    #--- convert back to original extent ---#
    rcl[idx] <- breaks_c$class
  }

  if (isTRUE(plot)) {
    data <- terra::as.data.frame(rastermetric)
    names(data) <- "val"

    nm <- as.character(names(rastermetric))

    data$var <- nm

    #--- plot histogram of mraster metric with associated break lines ---#

    p <- ggplot2::ggplot(data, ggplot2::aes(val)) +
      ggplot2::geom_histogram() +
      ggplot2::geom_vline(xintercept = breaks, linetype = "dashed") +
      ggplot2::ggtitle(paste0(nm, " histogram with defined breaks")) +
      ggplot2::xlab(paste0(nm))

    if (!is.null(mraster2)) {
      data2 <- terra::as.data.frame(rastermetric2)
      names(data2) <- "val"

      nm2 <- as.character(names(rastermetric2))

      data2$var <- nm2

      data2 <- rbind(data, data2)

      b1 <- data.frame(var = nm, brk = breaks)
      b2 <- data.frame(var = nm2, brk = breaks2)

      bs <- rbind(b1, b2)

      p <- ggplot2::ggplot(data2, ggplot2::aes(val)) +
        ggplot2::geom_histogram() +
        ggplot2::geom_vline(linetype = "dashed", data = bs, mapping = ggplot2::aes(xintercept = brk)) +
        ggplot2::facet_wrap(~var, scales = "free")
    }

    suppressMessages(print(p))

    #--- set colour palette ---#

    terra::plot(rcl, main = "User break defined strata", type = "classes")
  }

  #--- write file to disc ---#

  if (!is.null(filename)) {
    terra::writeRaster(rcl, filename, overwrite = overwrite, ...)
  }

  #--- Output based on 'details' to return raster alone or list with details ---#

  if (isTRUE(details)) {

    #--- output break points raster with associated breaks ---#
    breaks_rcl <- list(
      details = list(
        breaks = breaks,
        breaks2 = if (!missing(breaks2)) breaks2
      ),
      raster = rcl,
      plot = if (exists("p")) p
    )

    return(breaks_rcl)
  } else {

    #--- just output raster ---#

    return(rcl)
  }
}
