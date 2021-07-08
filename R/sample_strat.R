#' Stratified sampling
#' 
#' @description Sampling based on a stratified raster.
#' 
#' @family sample functions
#'
#' @inheritParams sample_srs
#' @param existing sf or data.frame.  Existing plot network.
#' @param include Logical. If \code{TRUE} include existing plots in \code{n} total.
#' @param wrow Numeric. Number of row in the focal window (default is 3).
#' @param wcol Numeric. Number of columns in the focal window (default is 3).
#' @param details Logical. If \code{FALSE} (default) output is sf object of 
#' stratified samples. If \code{TRUE} return a list
#' where \code{$details} additional sampling information and \code{$raster} 
#' is an sf object of stratified samples.
#' @param plot Logial. Plots existing (circles) and new (crosses) samples on the first band of mraster.
#'
#' @importFrom magrittr %>%
#' @importFrom methods is
#' 
#' @return An sf object with \code{n} stratified samples.
#' 
#' @export

sample_strat <- function(sraster,
                         n,
                         mindist = 100,
                         existing = NULL,
                         include = FALSE,
                         access = NULL,
                         buff_inner = NULL,
                         buff_outer = NULL,
                         wrow = 3,
                         wcol = 3,
                         plot = FALSE,
                         details = FALSE) 
{
  
  #--- Set global vars ---#
  x <- y <- NULL
  
  #--- Error management ---#
  if (!inherits(sraster, "SpatRaster"))
    stop("'sraster' must be type SpatRaster", call. = FALSE)
  
  if (any(! c("strata") %in% names(sraster)))
    stop("'sraster must have a layer named 'strata'")
  
  if (!is.numeric(mindist))
    stop("'mindist' must be type numeric")
  
  if (!is.numeric(n))
    stop("'n' must be type numeric")
  
  if (!is.numeric(wrow))
    stop("'wrow' must be type numeric")
  
  if (!is.numeric(wcol))
    stop("'wcol' must be type numeric")
  
  if (!is.logical(plot))
    stop("'plot' must be type logical")
  
  if (!is.logical(details))
    stop("'details' must be type logical")
  
  #--- if the sraster has multiple bands subset the band named strata ---#
  if(terra::nlyr(sraster) > 1){
    
    sraster <- terra::subset(sraster, "strata")
    
  }
  
  #--- determine crs of input sraster ---#
  crs <- terra::crs(sraster)
  
  #--- if existing samples are provided ensure they are in the proper format ---#
  
  if (is.null(existing)) {
    
    if (isTRUE(include))
      stop("'existing' must be provided when 'include' == TRUE")
    
    #--- if existing samples do not exist make an empty data.frame called addSamples ---#
    addSamples <- data.frame(strata = NA, X = NA, Y = NA)
    extraCols <- character(0)
    
  } else {
    
    #--- existing must be either a data.frame or an sf object with columns names 'X' 'Y' 'strata' ---#
    
    if (!inherits(existing, "data.frame") && !inherits(existing,"sf"))
      stop("'existing' must be a data.frame or sf object")
    
    if (any(! c("strata") %in% names(existing)) )
      stop("'existing' must have an attribute named 'strata'")
    
    #--- error handling in the presence of 'existing' ---#
    if (!inherits(existing,"sf"))
      stop("'existing' must be an 'sf' object")
    
      if(inherits(sf::st_geometry(existing),"sfc_POINT")){
        
        #--- if existing is an sf object extract the coordinates and the strata vector ---#
        
        exist_xy <- sf::st_coordinates(existing)
        
        strata <- existing$strata
        
        existing <- as.data.frame(cbind(strata, exist_xy))
      
      
    } else {
      
      stop("'existing' geometry type must be 'sfc_POINT'")
      
    }
    
    #--- if existing samples do exist ensure proper naming convention ---#
    
    if (any(! c("X", "Y") %in% colnames(existing)) ) {
      
      #--- if coordinate column names are lowercase change them to uppercase to match requirements ---#
      
      if (any(c("x", "y") %in% colnames(existing))) {
        
        existing <- existing %>%
          dplyr::rename(X = x,
                 Y = y)
        
        message("'existing' column coordinate names are lowercase - converting to uppercase")
        
      } else {
        
        #--- if no x/y columns are present stop ---#  
        
        stop("'existing' must have columns named 'X' and 'Y'")
        
      }
      
      
    }
    
    addSamples <- existing
    
  } 
  
  extraCols <- colnames(existing)[!colnames(existing) %in% c("X", "Y", "strata")]
  
  # Transform strata to numeric if factor
  if (is(addSamples$strata, "factor")) {
    addSamples$strata <- as.numeric(as.character(addSamples$strata))
    
  }
  
  #--- determine number of samples for each strata ---#
  
  if (isTRUE(include)) {
    message("'existing' samples being included in 'n' calculation")
    
    toSample <- tallySamples(sraster, n, existing)
    
  } else {
    
    toSample <- tallySamples(sraster,n)
    
  }
  
  
  #--- determine access buffers ---#
  
  if (!missing(access)){
    
    #--- error handling in the presence of 'access' ---#
    if (!inherits(access,"sf"))
      stop("'access' must be an 'sf' object")
      
    if(!inherits(sf::st_geometry(access),"sfc_MULTILINESTRING"))
     stop("'access' geometry type must be 'sfc_MULTILINESTRING'")
    
    if (buff_inner > buff_outer)
      stop("'buff_inner' must be < 'buff_outer'")
    
    #--- convert vectors to spatVector to synergize with terra sraster functions---#
    
    access <- terra::vect(access)
    
    #--- list all buffers to catch NULL values within error handling ---#
    buffers <- list(buff_inner, buff_outer)
    
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
    
    buff_in <- terra::buffer(x = access,
                             width = buff_inner)
    
    buff_out <- terra::buffer(x = access,
                              width = buff_outer)
    
    #--- make difference and aggregate inner and outer buffers to prevent sampling too close to access ---#
    
    buffer <- terra::aggregate(buff_out - buff_in)
    
    #--- mask sraster for plotting ---#
    raster_masked <- terra::mask(sraster, 
                                 mask = buffer)
    
  }
  
  ####################################
  #--- Start of sampling function ---#
  ####################################
  
  for (i in 1:nrow(toSample)) {
    s <- as.numeric(toSample[i, 1])
    nSamp <- as.numeric(toSample[i, 2])
    
    message(paste0("Processing strata : ", s))
    
    if (nSamp > 0) {
      #--- mask for individual strata ---#
      
      strata_m <- terra::mask(sraster,
                              mask = sraster,
                              maskvalues = s,
                              inverse = TRUE
      )
      names(strata_m) <- "strata"
      
      #--- if access line polygon is specified create inner and outer buffers
      
      if (!missing(access)) {
        
        strata_m_buff <- terra::mask(strata_m, 
                                     mask = buffer)
        
        sampAvail <- sum(!is.na(terra::values(strata_m_buff)))
        
        if (sampAvail > nSamp) {
          message(
            paste0(
              "Buffered area contains ",
              sampAvail,
              " available  candidates. Sampling to reach ",
              nSamp,
              " samples starting."
            )
          )
          
          #--- rename to original strata sraster that will be used for sampling ---#
          strata_m <- strata_m_buff
          
          #--- if there are no samples to take within the specified 'buff_outer' distance extend buffer until values are found ---#
          
        } else {
          
          stop("Insufficient candidate samples within the buffered access extent. Consider altering buffer widths.")
          
        }
        
      }
      
      ##################
      #--- sampling ---#
      ##################
      
      #--- RULE 1: select only cells surrounded by cells with same strata ---#
      
      #--- Define focal window ---#
      w <- matrix(1 / (wrow * wcol), nr = wrow, nc = wcol)
      
      suppressWarnings(strata_m_clust <-
                         terra::focal(
                           strata_m,
                           w = w,
                           na.rm = FALSE,
                           na.only = FALSE
                         ))
      names(strata_m_clust) <- "strata"
      
      #--- Initiate number of sampled cells ---#
      add_strata <- addSamples %>%
        dplyr::filter(strata == s)
      
      if (nrow(add_strata) > 0) {
        add_strata$type <- "existing"
        
        if (!"rule" %in% colnames(add_strata)) {
          add_strata$rule <- "existing"
          
        }
      }
      
      #--- create indices for all, NA, and valid sampling candidates ---#
      
      idx_all <- 1:terra::ncell(strata_m_clust)
      idx_na <- is.na(terra::values(strata_m_clust))
      validCandidates <- idx_all[!idx_na]
      
      #--- Rule 1 sampling ---#
      nCount <- 0 #Number of sampled cells
      
      # While loop for RULE 1
      while (length(validCandidates) > 0 & nCount < nSamp) {
        #-- identify potential sample from candidates ---#
        smp <- sample(1:length(validCandidates), size = 1)
        
        smp_cell <- validCandidates[smp]
        
        #--- extract coordinates and sample details ---#
        
        add_temp <- data.frame(
          cell = smp_cell,
          X = terra::xFromCell(strata_m_clust, smp_cell),
          Y = terra::yFromCell(strata_m_clust, smp_cell),
          strata = strata_m_clust[smp_cell]
        )
        
        #--- Remove sampled cell from validCandidates so that it cannot be sampled again later ---#
        
        validCandidates <- validCandidates[-smp]
        
        #--- populate add_temp with values ---#
        add_temp$type <- "new"
        add_temp$rule <- "rule1"
        add_temp[, extraCols] <- NA
        
        #--- If add_strata is empty, sampled cell accepted ---#
        
        if (nrow(add_strata) == 0) {
          
          add_strata <- add_temp[, c("X", "Y", "strata", "type", "rule", extraCols)]
          
          nCount <-  nCount + 1
          
          #--- If add_strata isn't empty, check distance with all other sampled cells in strata ---#
        } else {
          
          dist <- spatstat.geom::crossdist(add_temp$X, add_temp$Y , add_strata$X , add_strata$Y)
          
          #--- If all less than 'mindist' - accept sampled cell otherwise reject ---#
          if (all(as.numeric(dist) > mindist)) {
            
            add_strata <- rbind(add_strata, add_temp[, c("X", "Y", "strata", "type", "rule", extraCols)])
            
            nCount <-  nCount + 1
            
          }
        }
      }
      
      #---- RULE 3 sampling ---#
      
      if (nCount < nSamp) {
        
        idx_all <- 1:terra::ncell(strata_m)
        idx_na <- is.na(terra::values(strata_m))
        validCandidates <- idx_all[!idx_na]
        
        while(length(validCandidates) > 0 & nCount < nSamp){
          
          #-- identify potential sample from candidates ---#
          smp <- sample(1:length(validCandidates), size = 1)
          
          smp_cell <- validCandidates[smp]
          
          #--- extract coordinates and sample details ---#
          
          add_temp <- data.frame(
            cell = smp_cell,
            X = terra::xFromCell(strata_m, smp_cell),
            Y = terra::yFromCell(strata_m, smp_cell),
            strata = validCandidates[smp_cell]
          )
          
          #--- Remove sampled cell from validCandidates so that it cannot be sampled again later ---#
          
          validCandidates <- validCandidates[-smp]
          
          add_temp$rule <- "isolated"
          add_temp$type <- "new"
          add_temp[,extraCols] <- NA
          add_temp$strata <- s
          
          if (nrow(add_strata) == 0) {
            
            add_strata <- add_temp[,c("X","Y","strata","type","rule",extraCols)]
            
            nCount <- nCount + 1
            
          } else {
            
            dist <- spatstat.geom::crossdist(add_temp$X,add_temp$Y,add_strata$X,add_strata$Y)
            
            if( all(as.numeric(dist) > mindist )){
              
              add_strata <- rbind(add_strata, add_temp[,c("X","Y","strata","type","rule",extraCols)])
              
              nCount <-  nCount + 1
              
            }
          }
        }
      }
      
      if (nCount < nSamp){
        
        message(sprintf("Strata %s: couldn't select required number of samples: %i instead of %i \n", s , nCount , nSamp))
        
      }
      
    }
    
    # Create out object if first iteration of loop
    # Else just rbind output with what has been processed in the loop
    if (i == 1) {
      
      out <- add_strata
      
    }else{
      
      out <- rbind(out,add_strata)
      
    }
    
  }
  
  #--- convert coordinates to a spatial points object ---#
  samples <- out %>%
    as.data.frame() %>%
    sf::st_as_sf(., coords = c("X", "Y"))
  
  #--- assign sraster crs to spatial points object ---#
  sf::st_crs(samples) <- crs
  
  
  #--- plot the raster and samples if desired ---#
  
  if(isTRUE(plot)){
    
    #--- if existing is not provided plot the masked raster ---#
    
    if(missing(existing)){
      
      #--- if access is also missing plot the full sraster extent ---#
      
      if(missing(access)){
        
        terra::plot(sraster)
        suppressWarnings(terra::plot(samples, add = T, col = "black", pch = ifelse(samples$type=="existing",1,3)))
        
        #--- if access is provided plot the masked access sraster ---#
        
      } else {

        terra::plot(raster_masked)
        suppressWarnings(terra::plot(samples, add = T, col = "black", pch = ifelse(samples$type=="existing",1,3)))
      
      }
      
      #--- if existing is provided plot the full raster ---#
      
    } else {
    
      #--- plot input sraster and random samples ---#
      
      terra::plot(sraster)
      suppressWarnings(terra::plot(samples, add = T, col = "black", pch = ifelse(samples$type=="existing",1,3)))
    
    }
    

  }

  if ( isTRUE(details) ){
    
    #--- output metrics details along with stratification raster ---#
    
    output <- list(sampleDist = toSample, samples = samples)
    
    #--- output samples dataframe ---#
    return(output)
    
    
  } else {
    
    #--- just output raster ---#
    
    return(samples)
    
  
    }
  
}
