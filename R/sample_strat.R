# raster = spatRaster. Single band stratification raster derived from strat_* functions.
# ns = Numeric. Number of desired samples.
# mindist = Numeric. Minimum allowable distance between selected samples.
# existing = data.frame or sf. Existing sample plot network. must have columns denoted X and Y.
# access = spatVector. Road access network - must be lines.
# buff_inner = Numeric. Inner buffer boundary specifying distance from access where plots cannot be sampled.
# buff_outer = Numeric. Outer buffer boundary specifying distance from access where plots can be sampled.
# buff_extend = Numeric. Amount 'buff_outer' will be increased by if no valid sample candidates are available.
# buff_max = Numeric. Maximum buffer distance resulting from increasing iterations of 'buff_outer' + 'buff_extend'.
# wrow = Numeric. Number of row in the focal window (default is 3).
# wcol = Numeric. Number of columns in the focal window (default is 3).

sample_strat <- function(raster,
                         ns,
                         mindist,
                         existing = NULL,
                         access = NULL,
                         buff_inner = NULL,
                         buff_outer = NULL,
                         buff_extend = NULL,
                         buff_max = NULL,
                         wrow = 3,
                         wcol = 3) {
  #--- Error management ---#
  if (!inherits(raster, "SpatRaster"))
    stop("'raster' must be type SpatRaster", call. = FALSE)
  
  if (!is.numeric(mindist))
    stop("'mindist' must be type numeric")
  
  if (!is.numeric(ns))
    stop("'ns' must be type numeric")
  
  if (!is.numeric(wrow))
    stop("'wrow' must be type numeric")
  
  if (!is.numeric(wcol))
    stop("'wcol' must be type numeric")
  
  #--- determine crs of input raster ---#
  crs <- crs(raster)
  
  #--- if existing samples are provided ensure they are in the proper format ---#
  
  if (is.null(existing)) {
    #--- if existing samples do not exist make an empty data.frame called addSamples ---#
    addSamples <- data.frame(strata = NA, X = NA, Y = NA)
    extraCols <- character(0)
    
  } else {
    
    #--- existing must be either a data.frame or an sf object with columns names 'X' 'Y' 'strata' ---#
    
    if (!inherits(existing, "data.frame") && !inherits(existing,"sf"))
      stop("'existing' must be a data.frame or sf object")
    
    if (any(! c("strata") %in% names(existing)) )
      stop("'existing' must have an attribute named 'strata'")
    
    if (inherits(existing,"sf") && inherits(sf::st_geometry(existing),"sfc_POINT")){
      
      #--- if existing is an sf object extract the coordinates and the strata vector ---#
      
      exist_xy <- st_coordinates(existing)
      
      strata <- existing$strata
      
      existing <- as.data.frame(cbind(strata, exist_xy))

      
    }
    
    #--- if existing samples do exist ensure proper naming convention ---#
    
    if (any(! c("X", "Y") %in% colnames(existing)) ) {
      
      #--- if coordinate column names are lowercase change them to uppercase to match requirements ---#
      
      if (any(c("x", "y") %in% colnames(existing))) {
        
        existing <- existing %>%
          rename(X = x,
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
  
  toSample <- tallySamples(raster, ns)
  
  #--- determine access buffers ---#
  
  for (i in 1:nrow(toSample)) {
    s <- as.numeric(toSample[i, 1])
    nsamp <- as.numeric(toSample[i, 2])
    
    print(paste0("Processing strata : ", s))
    
    if (nsamp > 0) {
      #--- mask for individual strata ---#
      
      strata_m <- terra::mask(raster,
                             mask = raster,
                             maskvalues = s,
                             inverse = TRUE
      )
      names(strata_m) <- "strata"
      
      #--- if access line polygon is specified create inner and outer buffers
      
      if (!is.null(access)) {
        
        #--- error handling in the presence of 'access' ---#
        if (!inherits(access, "SpatVector"))
          stop("'access' must be type SpatVector", call. = FALSE)
        
        if (any(buff_max < c(buff_inner, buff_outer)))
          stop("'buff_inner' must be < 'buff_outer' & 'buff_outer' must be < 'buff_max'")

        #--- list all buffers to catch NULL values within error handling ---#
        buffers <- list(buff_inner, buff_outer, buff_extend, buff_max)

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
        
        strata_m_buff <- terra::mask(strata_m, 
                                     mask = buffer)
        
        sampAvail <- sum(!is.na(values(strata_m_buff)))
        
        if (sampAvail > nsamp) {
          message(
            paste0(
              "Buffered area contains ",
              sampAvail,
              " available  candidates. Sampling to reach ",
              nsamp,
              " samples starting."
            )
          )
          
          #--- rename to original strata raster that will be used for sampling ---#
          strata_m <- strata_m_buff
          
          #--- if there are no samples to take within the specified 'buff_outer' distance extend buffer until values are found ---#
          
        } else {
          if (is.null(buff_extend))
            stop("Insufficient candidate samples are within the buffer extent and 'buff_extend' is null. Supply a value to iterate increasing external buffers or increase 'buff_outer' value."
            )
          
          counter <- 1
          
          while (sampAvail < nsamp) {
            #--- extend buffer based on 'buff_extend' ---#
            buff_outer_n <- buff_outer + (buff_extend * counter)
            
            #--- if the max buffer has been reached stop ---#
            if (buff_outer_n > buff_max) {
              stop(
                "Maximum buffer size of 'buff_max' has been reached. Insufficient number of candidates to reach sample size."
              )
              
            } else {
              message(
                paste0(
                  "Buffered area contains no available samples. Increasing buff_outer to ",
                  buff_outer_n,
                  " m"
                )
              )
              
            }
            
            #--- recompute outer buffer with buffer extension ---#
            buff_out <-
              terra::buffer(x = roads,
                            width = buff_outer_n,
                            capstyle = "round")
            
            buffer <- aggregate(buff_out - buff_in)
            
            #--- mask strata raster with extended buffer ---#
            strata_m_buff <- terra::mask(strata_m, mask = buffer)
            
            sampAvail <- sum(!is.na(values(strata_m_buff)))
            
            #--- if number of  candidate samples > samples needed exit while loop and begin sampling ---#
            if (sampAvail > nsamp) {
              message(
                paste0(
                  "External buffer of ",
                  buff_outer_n,
                  " m contains ",
                  sampAvail,
                  " available  candidates. Sampling to reach ",
                  nsamp,
                  " samples starting."
                )
              )
              
              #--- rename to original strata raster that will be used for sampling ---#
              strata_m <- strata_m_buff
              
            }
            counter <- counter + 1
            
          }
          
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
        add_strata$type <- "Existing"
        
        if (!"rule" %in% colnames(add_strata)) {
          add_strata$rule <- NA
          
        }
      }
      
      #--- create indices for all, NA, and valid sampling candidates ---#
      
      idx_all <- 1:ncell(strata_m_clust)
      idx_na <- is.na(terra::values(strata_m_clust))
      validCandidates <- idx_all[!idx_na]
      
      #--- Rule 1 sampling ---#
      nCount <- 0 #Number of sampled cells
      
      # While loop for RULE 1
      while (length(validCandidates) > 0 & nCount < nsamp) {
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
        add_temp$type <- "New"
        add_temp$rule <- "rule1"
        add_temp[, extraCols] <- NA
        
        #--- If add_strata is empty, sampled cell accepted ---#
        
        if (nrow(add_strata) == 0) {
          
          add_strata <- add_temp[, c("X", "Y", "strata", "type", "rule", extraCols)]
          
          nCount <-  nCount + 1
          
          #--- If add_strata isn't empty, check distance with all other sampled cells in strata ---#
        } else {
          
          dist <- crossdist(add_temp$X, add_temp$Y , add_strata$X , add_strata$Y)
          
          #--- If all less than 'mindist' - accept sampled cell otherwise reject ---#
          if (all(as.numeric(dist) > mindist)) {
            
            add_strata <- rbind(add_strata, add_temp[, c("X", "Y", "strata", "type", "rule", extraCols)])
            
            nCount <-  nCount + 1
            
          }
        }
      }
      
      #---- RULE 3 sampling ---#
      
      if (nCount < nsamp) {
        
        idx_all <- 1:ncell(strata_m)
        idx_na <- is.na(terra::values(strata_m))
        validCandidates <- idx_all[!idx_na]
        
        while(length(validCandidates) > 0 & nCount < nsamp){
          
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
          add_temp$type <- "New"
          add_temp[,extraCols] <- NA
          add_temp$strata <- s
          
          if (nrow(add_strata) == 0) {
            
            add_strata <- add_temp[,c("X","Y","strata","type","rule",extraCols)]
            
            nCount = nCount +1
            
          } else {
            
            dist <- crossdist(add_temp$X,add_temp$Y,add_strata$X,add_strata$Y)
            
            if( all(as.numeric(dist) > mindist )){
              
              add_strata <- rbind(add_strata, add_temp[,c("X","Y","strata","type","rule",extraCols)])
              
              nCount <-  nCount + 1
              
            }
          }
        }
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
    st_as_sf(., coords = c("X", "Y"))
  
  #--- assign raster crs to spatial points object ---#
  st_crs(samples) <- crs
  
  #--- plot input raster and random samples ---#
  terra::plot(raster[[1]])
  suppressWarnings(terra::plot(samples, add = T, col = ifelse(samples$type=="Existing","Red","Black")),)
  
  #--- output samples dataframe ---#
  return(samples)
  
}
