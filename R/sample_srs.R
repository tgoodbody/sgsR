#' Simple random sampling
#'
#' @description Randomly sample within a stratification raster extent.
#'
#' @family sample functions
#'
#' @inheritParams sample_systematic
#'
#' @param raster spatRaster. Raster to be used for random sampling.
#' @param nSamp Numeric. Number of desired samples.
#' @param mindist Numeric. Minimum allowable distance between selected
#'  samples. \code{Default = NULL}.
#' @param access sf 'LINESTRING' or 'MULTILINESTRING'. Access network.
#' @param buff_inner Numeric. Inner buffer boundary specifying distance
#'  from access where plots cannot be sampled.
#' @param buff_outer Numeric. Outer buffer boundary specifying distance
#'  from access where plots can be sampled.
#' @param plot Logical. Plots output strata raster with samples.
#'
#' @return An sf object with \code{nSamp} randomly sampled points.
#'
#' @examples
#' #--- Load raster and access files ---#
#' r <- system.file("extdata", "sraster.tif", package = "sgsR")
#' sr <- terra::rast(r)
#'
#' a <- system.file("extdata", "access.shp", package = "sgsR")
#' ac <- sf::st_read(a)
#'
#' #--- perform simple random sampling ---#
#' sample_srs(
#'   raster = sr,
#'   nSamp = 200,
#' )
#'
#' @author Tristan R.H. Goodbody & Martin Queinnec
#'
#' @export

sample_srs <- function(raster,
                       nSamp,
                       mindist = NULL,
                       access = NULL,
                       buff_inner = NULL,
                       buff_outer = NULL,
                       plot = FALSE,
                       filename = NULL,
                       overwrite = FALSE) {
  #--- Error management ---#

  if (!inherits(raster, "SpatRaster")) {
    stop("'raster' must be type SpatRaster.", call. = FALSE)
  }

  if (!is.numeric(nSamp)) {
    stop("'nSamp' must be type numeric.", call. = FALSE)
  }

  if (!is.null(mindist)) {
    if (!is.numeric(mindist)) {
      stop("'mindist' must be type numeric.", call. = FALSE)
    }
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical.", call. = FALSE)
  }

  if (is.na(terra::crs(raster))) {
    stop("'raster' does not have a coordinate system.", call. = FALSE)
  }

  rasterP <- raster <- raster[[1]]

  #--- determine crs of input raster ---#
  crs <- terra::crs(raster)

  if (!is.null(access)) {
    access_buff <- mask_access(raster = raster, access = access, buff_inner = buff_inner, buff_outer = buff_outer)

    raster <- access_buff$rast
  }

  #--- create empty dataframe for samples to be populated to ---#
  add_strata <- data.frame()

  #--- create indices for all, NA, and valid sampling candidates ---#

  idx_all <- 1:terra::ncell(raster)
  idx_na <- is.na(terra::values(raster))
  validCandidates <- idx_all[!idx_na]

  #--- Rule 1 sampling ---#
  nCount <- 0 # Number of sampled cells
  iter <- 0 # Number of sample iterations

  # While loop for RULE 1
  while (length(validCandidates) > 0 & nCount < nSamp) {
    #-- identify potential sample from candidates ---#
    smp <- sample(1:length(validCandidates), size = 1)

    smp_cell <- validCandidates[smp]

    #--- extract coordinates and sample details ---#

    add_temp <- data.frame(
      cell = smp_cell,
      X = terra::xFromCell(raster, smp_cell),
      Y = terra::yFromCell(raster, smp_cell)
    )

    #--- Remove sampled cell from validCandidates so that it cannot be sampled again later ---#

    validCandidates <- validCandidates[-smp]

    #--- If add_strata is empty, sampled cell accepted ---#

    if (nrow(add_strata) == 0) {
      add_strata <- add_temp[, c("X", "Y")]

      nCount <- nCount + 1

      #--- If add_strata isnt empty, check distance with all other sampled cells in strata ---#
    }

    if (!is.null(mindist)) {
      dist <- spatstat.geom::crossdist(add_temp$X, add_temp$Y, add_strata$X, add_strata$Y)

      #--- If all less than 'mindist' - accept sampled cell otherwise reject ---#
      if (all(as.numeric(dist) > mindist)) {
        add_strata <- rbind(add_strata, add_temp[, c("X", "Y")])

        nCount <- nCount + 1
      }
    } else if (iter != 0) {
      add_strata <- rbind(add_strata, add_temp[, c("X", "Y")])

      nCount <- nCount + 1
    }

    iter <- iter + 1
  }

  if (nrow(add_strata) < nSamp) {
    message(paste0("Sampling was not able to select ", nSamp, " sample units. Output has ", nrow(add_strata), " sample units."))
  }

  #--- convert coordinates to a spatial points object ---#
  samples <- add_strata %>%
    as.data.frame() %>%
    sf::st_as_sf(., coords = c("X", "Y"), crs = crs)

  if (isTRUE(plot)) {
    if (!is.null(access)) {
      #--- plot input raster and random samples ---#
      terra::plot(rasterP[[1]])
      suppressWarnings(terra::plot(access_buff$buff, add = T, border = c("gray30"), col = "gray10", alpha = 0.1))
      suppressWarnings(terra::plot(samples, add = T, col = "black"))
    } else {
      #--- plot input raster and random samples ---#
      terra::plot(rasterP[[1]])
      suppressWarnings(terra::plot(samples, add = T, col = "black"))
    }
  }

  #--- write outputs if desired ---#
  write_samples(samples = samples, filename = filename, overwrite = overwrite)

  #--- output samples sf ---#

  return(samples)
}
