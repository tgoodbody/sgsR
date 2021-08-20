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
#'  samples. Default = 100.
#' @param access sf. Road access network - must be lines.
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
#' r <- system.file("extdata", "kmeans.tif", package = "sgsR")
#' sr <- terra::rast(r)
#'
#' a <- system.file("extdata", "roads.shp", package = "sgsR")
#' ac <- sf::st_read(a)
#'
#' #--- perform simple random sampling ---#
#' sample_srs(raster = sr, 
#'            nSamp = 200, 
#'            plot = TRUE)
#' 
#' sample_srs(raster = sr, 
#'            nSamp = 200, 
#'            access = ac,
#'            mindist = 200,
#'            buff_inner = 50,
#'            buff_outer = 200)
#' 
#' sample_srs(raster = sr,
#'            nSamp = 200,
#'            access = ac,
#'            buff_inner = 50,
#'            buff_outer = 200,
#'            filename = tempfile(fileext = ".shp"))
#'            
#' @author Tristan R.H. Goodbody & Martin Queinnec
#'
#' @export

sample_srs <- function(raster,
                       nSamp,
                       mindist = 100,
                       access = NULL,
                       buff_inner = NULL,
                       buff_outer = NULL,
                       plot = FALSE,
                       filename = NULL,
                       overwrite = FALSE) {

  #--- Error management ---#

  if (!inherits(raster, "SpatRaster")) {
    stop("'raster' must be type SpatRaster", call. = FALSE)
  }

  if (!is.numeric(nSamp)) {
    stop("'nSamp' must be type numeric")
  }

  if (!is.numeric(mindist)) {
    stop("'mindist' must be type numeric")
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical")
  }

  ######################################
  ## DETERMINE NULL / NA SYNTAX FOR CRS##
  ######################################

  if (is.na(terra::crs(raster))) {
    stop("'raster' does not have a coordinate system")
  }

  rasterP <- raster <- raster[[1]]

  #--- determine crs of input raster ---#
  crs <- terra::crs(raster, proj = TRUE)

  if (!is.null(access)) {

    #--- error handling in the presence of 'access' ---#
    if (!inherits(access, "sf")) {
      stop("'access' must be an 'sf' object")
    }

    if (!inherits(sf::st_geometry(access), "sfc_MULTILINESTRING")) {
      stop("'access' geometry type must be 'sfc_MULTILINESTRING'")
    }

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

      #--- If awdd_strata isnt empty, check distance with all other sampled cells in strata ---#
    }

    if (!is.null(mindist)) {
      dist <- spatstat.geom::crossdist(add_temp$X, add_temp$Y, add_strata$X, add_strata$Y)

      #--- If all less than 'mindist' - accept sampled cell otherwise reject ---#
      if (all(as.numeric(dist) > mindist)) {
        add_strata <- rbind(add_strata, add_temp[, c("X", "Y")])

        nCount <- nCount + 1
      }
    } else {
      add_strata <- rbind(add_strata, add_temp[, c("X", "Y")])

      nCount <- nCount + 1
    }
  }

  #--- convert coordinates to a spatial points object ---#
  samples <- add_strata %>%
    as.data.frame() %>%
    sf::st_as_sf(., coords = c("X", "Y"))

  #--- assign raster crs to spatial points object ---#
  sf::st_crs(samples) <- crs

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

  if (!is.null(filename)) {
    if (!is.logical(overwrite)) {
      stop("'overwrite' must be either TRUE or FALSE")
    }

    if (file.exists(filename) & isFALSE(overwrite)) {
      stop(paste0(filename, " already exists and overwrite = FALSE"))
    }

    sf::st_write(samples, filename, delete_layer = overwrite)
  }

  #--- output samples sf ---#

  return(samples)
}
