# raster = spatRaster. Raster to be used for simple random sampling
# ns = Numeric. Number of desired samples.
# mindist = Numeric. Minimum allowable distance between selected samples.
# access = spatVector. Road access network - must be lines.
# buff_inner = Numeric. Inner buffer boundary specifying distance from access where plots cannot be sampled.
# buff_outer = Numeric. Outer buffer boundary specifying distance from access where plots can be sampled.


sample_srs <- function(raster,
                       ns,
                       mindist,
                       access = NULL,
                       buff_inner = NULL,
                       buff_outer = NULL) {
  #--- Error management ---#
  if (!inherits(raster, "SpatRaster"))
    stop("'raster' must be type SpatRaster", call. = FALSE)

  if (!is.numeric(ns))
    stop("'ns' must be type numeric")

  if (!is.numeric(mindist))
    stop("'mindist' must be type numeric")
  
  ######################################
  ##DETERMINE NULL / NA SYNTAX FOR CRS##
  ######################################

  if (is.na(crs(raster)))
    stop("'raster' does not have a coordinate system")

  raster <- raster[[1]]
  
  #--- determine crs of input raster ---#
  crs <- crs(raster)

  if (!is.null(access)) {

    if (!inherits(access, "SpatVector"))
      stop("'access' must be type SpatVector", call. = FALSE)

    #--- list all buffers to catch NULL values within error handling ---#
    buffers <- list(buff_inner, buff_outer)

    #--- error handling in the presence of 'access' ---#
    if (any(vapply(buffers, is.null, TRUE)))
      stop("All 'buff_*' paramaters must be provided when 'access' is defined.")

    if (!any(vapply(buffers, is.numeric, FALSE)))
      stop("All 'buff_*' paramaters must be type numeric")

    message(
      paste0(
        "An access layer has been provided. An internal buffer of ",
        buff_inner,
        " m and an external buffer of ",
        buff_outer,
        " m have been applied"
      )
    )

    #--- make access buffer with user defined values ---#

    buff_in <- terra::buffer(x = roads,
                             width = buff_inner,
                             capstyle = "round")

    buff_out <- terra::buffer(x = roads,
                              width = buff_outer,
                              capstyle = "round")

    #--- make difference and aggregate inner and outer buffers to prevent sampling too close to access ---#
    buffer <- aggregate(buff_out - buff_in)

    raster <- terra::mask(raster, mask = buffer)

  }
  
  #--- create empty dataframe for samples to be populated to ---#
  add_strata <- data.frame()
  
  #--- create indices for all, NA, and valid sampling candidates ---#

  idx_all <- 1:ncell(raster)
  idx_na <- is.na(terra::values(raster))
  validCandidates <- idx_all[!idx_na]

  #--- Rule 1 sampling ---#
  nCount <- 0 #Number of sampled cells

  # While loop for RULE 1
  while (length(validCandidates) > 0 & nCount < ns) {
    #-- identify potential sample from candidates ---#
    smp <- sample(1:length(validCandidates), size = 1)

    smp_cell <- validCandidates[smp]

    #--- extract coordinates and sample details ---#

    add_temp <- data.frame(
      cell = smp_cell,
      x = terra::xFromCell(raster, smp_cell),
      y = terra::yFromCell(raster, smp_cell)
    )

    #--- Remove sampled cell from validCandidates so that it cannot be sampled again later ---#

    validCandidates <- validCandidates[-smp]

    #--- If add_strata is empty, sampled cell accepted ---#

    if (nrow(add_strata) == 0) {

      add_strata <- add_temp[, c("x", "y")]

      nCount <-  nCount + 1

      #--- If add_strata isnt empty, check distance with all other sampled cells in strata ---#
    } else {

      dist <- crossdist(add_temp$x, add_temp$y , add_strata$x , add_strata$y)

      #--- If all less than 'mindist' - accept sampled cell otherwise reject ---#
      if (all(as.numeric(dist) > mindist)) {

        add_strata <- rbind(add_strata, add_temp[, c("x", "y")])

        nCount <-  nCount + 1

      }
    }
  }
  
    #--- convert coordinates to a spatial points object ---#
    samples <- add_strata %>%
      as.data.frame() %>%
      st_as_sf(., coords = c("x", "y"))

    #--- assign raster crs to spatial points object ---#
    st_crs(samples) <- crs

    #--- plot input raster and random samples ---#
    terra::plot(raster[[1]])
    suppressWarnings(terra::plot(samples, add = T, col = "black"))

    #--- output samples dataframe ---#
    return(samples)


}


