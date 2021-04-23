#' Balanced raster sampling using \code{\link{BalancedSampling}} and \code{\link{SamplingBigData}} methods
#' @family sample functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams sample_srs
#' @param algorithm Character. One of \code{lpm2 lcube lcubestratified}
#' @param p Numeric. Inclusion probability for each candidate sample. Default is \code{n / N}
#' 
#' 
#' @return An sf object with \code{n} randomly sampled points.
#' 
#' @export

sample_balanced <- function(mraster,
                            n,
                            algorithm = "lpm2",
                            p = NULL,
                            access = NULL,
                            buff_inner = NULL,
                            buff_outer = NULL,
                            plot = FALSE) {
  
  #--- Error management ---#
  if (!inherits(mraster, "SpatRaster"))
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  
  if (!is.numeric(n))
    stop("'n' must be type numeric")
  
  if (!is.logical(plot))
    stop("'plot' must be type logical")
  
  if (!is.character(algorithm))
    stop("'algorith' must be type character")
  
  #--- list all available algorithms to determine if a valid one has been supplied ---#
  algs <- c("lpm2", "lcube", "lcubestratified")
  
  if (!algorithm %in% algs)
    stop("Unknown algorithm specified. Please use one of 'lpm2' 'lcube' 'lcubestratified'")
  
  ######################################
  ##DETERMINE NULL / NA SYNTAX FOR CRS##
  ######################################
  
  if (is.na(crs(mraster)))
    stop("'mraster' does not have a coordinate system")
  
  #--- determine crs of input mraster ---#
  crs <- crs(mraster)
  
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
    roads <- terra::vect(roads)
    
    #--- make access buffer with user defined values ---#
    
    buff_in <- terra::buffer(x = roads,
                             width = buff_inner,
                             capstyle = "round")
    
    buff_out <- terra::buffer(x = roads,
                              width = buff_outer,
                              capstyle = "round")
    
    #--- make difference and aggregate inner and outer buffers to prevent sampling too close to access ---#
    buffer <- terra::aggregate(buff_out - buff_in)
    
    mraster <- terra::mask(mraster, mask = buffer)
    
  }
  
  #--- extract XY coordinates from raster ---#
  vals <- terra::as.data.frame(mraster,xy=TRUE) %>%
    rename(X = x,
           Y = y)
  
  #--- # drop x,y matrix of auxiliary variables ---#
  vals_m = as.matrix(dplyr::select(vals,-X, -Y)) 
  
  N <- nrow(vals)
  
  #--- inclusion probability ---#
  
  if (is.null(p)) {
    
    #--- if 'p' is not defined use the default ---#
    
    p <- rep(n/N,N)
    
  } else {
    
    if (!is.numeric(p))
      stop("'p' must be type numeric")
    
  }
  
  
  if (algorithm == "lpm2"){
    
    sampled <- SamplingBigData::lpm2_kdtree(prob = p, x = vals_m)
    
    
   }
  
  if (algorithm == "lcube"){
    
    sampled <- BalancedSampling::lcube(prob = p, Xspread = vals_m, Xbal = cbind(p))
    
    
  }
  
  if (algorithm == "lcubestratified"){
    
    if(!"strata" %in% names(mraster))
      stop("'mraster' must have a variable named 'strata' to use the 'lcubestratified' algorithm")
    
    #--- create indices for all, NA, and valid sampling candidates ---#
    
    strata_v <- as.vector(vals$strata)
    
    #--- remove strata as a sampling variable and convert to matrix ---#
    
    #--- # drop x,y matrix of auxiliary variables ---#
    vals_m = as.matrix(dplyr::select(vals,-X, -Y, -strata)) 

    #--- generates a binary output where 0 is not sampled and 1 is sampled ---#
    
    sampled <- BalancedSampling::lcubestratified(prob = p, 
                                                 Xspread = vals_m, 
                                                 Xbal = cbind(p),
                                                 integerStrata = strata_v)
    
    #--- extract all 1 (sampled) cells ---#
    
    sampled <-  (1:N)[sampled==1]
    
    
  }
  
  samples <- vals[sampled,]

  #--- convert coordinates to a spatial points object ---#
  samples <- dplyr::select(samples, X, Y) %>%
    as.data.frame() %>%
    st_as_sf(., coords = c("X", "Y"))
  
  #--- assign mraster crs to spatial points object ---#
  st_crs(samples) <- crs
  
  if(isTRUE(plot)){
    
    #--- plot input mraster and random samples ---#
    terra::plot(mraster[[1]])
    suppressWarnings(terra::plot(samples, add = T, col = "black"))
    
  }
  
  #--- output samples sf ---#
  
  return(samples)
  
}



