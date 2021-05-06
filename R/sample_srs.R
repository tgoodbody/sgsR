#' Randomly sample within a stratification raster extent.
#' @family sample functions
#'
#' @param sraster spatRaster. Stratification raster to be used for sampling.
#' @param n Numeric. Number of desired samples.
#' @param mindist Numeric. Minimum allowable distance between selected
#'  samples. Default = 100.
#' @param access sf. Road access network - must be lines.
#' @param buff_inner Numeric. Inner buffer boundary specifying distance
#'  from access where plots cannot be sampled.
#' @param buff_outer Numeric. Outer buffer boundary specifying distance
#'  from access where plots can be sampled.
#' @param plot Logical. Plots output strata raster with samples.
#' 
#' @return An sf object with \code{n} randomly sampled points.
#' 
#' @export

sample_srs <- function(sraster,
                       n,
                       mindist = 100,
                       access = NULL,
                       buff_inner = NULL,
                       buff_outer = NULL,
                       plot = FALSE) {
  #--- Error management ---#
  if (!inherits(sraster, "SpatRaster"))
    stop("'sraster' must be type SpatRaster", call. = FALSE)

  if (!is.numeric(n))
    stop("'n' must be type numeric")

  if (!is.numeric(mindist))
    stop("'mindist' must be type numeric")
  
  if (!is.logical(plot))
    stop("'plot' must be type logical")
  
  ######################################
  ##DETERMINE NULL / NA SYNTAX FOR CRS##
  ######################################

  if (is.na(crs(sraster)))
    stop("'sraster' does not have a coordinate system")

  sraster <- sraster[[1]]
  
  #--- determine crs of input sraster ---#
  crs <- crs(sraster)

  if (!is.null(access)) {

    if (!inherits(access,"sf") && inherits(sf::st_geometry(access),"sfc_MULTILINESTRING"))
      stop("'access' must be an 'sf' object of type 'sfc_MULTILINESTRING' geometry")
    
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

    #--- convert vectors to spatVector to synergize with terra raster functions---#
    roads <- vect(roads)
    
    #--- make access buffer with user defined values ---#

    buff_in <- terra::buffer(x = roads,
                             width = buff_inner)

    buff_out <- terra::buffer(x = roads,
                              width = buff_outer)

    #--- make difference and aggregate inner and outer buffers to prevent sampling too close to access ---#
    buffer <- aggregate(buff_out - buff_in)

    sraster <- terra::mask(sraster, mask = buffer)

  }
  
  #--- create empty dataframe for samples to be populated to ---#
  add_strata <- data.frame()
  
  #--- create indices for all, NA, and valid sampling candidates ---#

  idx_all <- 1:ncell(sraster)
  idx_na <- is.na(terra::values(sraster))
  validCandidates <- idx_all[!idx_na]

  #--- Rule 1 sampling ---#
  nCount <- 0 #Number of sampled cells

  # While loop for RULE 1
  while (length(validCandidates) > 0 & nCount < n) {
    #-- identify potential sample from candidates ---#
    smp <- sample(1:length(validCandidates), size = 1)

    smp_cell <- validCandidates[smp]

    #--- extract coordinates and sample details ---#

    add_temp <- data.frame(
      cell = smp_cell,
      x = terra::xFromCell(sraster, smp_cell),
      y = terra::yFromCell(sraster, smp_cell)
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

    #--- assign sraster crs to spatial points object ---#
    st_crs(samples) <- crs
    
    if(isTRUE(plot)){

    #--- plot input sraster and random samples ---#
    terra::plot(sraster[[1]])
    suppressWarnings(terra::plot(samples, add = T, col = "black"))

    }
      
      #--- output samples sf ---#
      
      return(samples)

}


