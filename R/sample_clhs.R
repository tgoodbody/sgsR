#' Conditioned Latin Hypercube Sampling using the \code{clhs} package functionality
#' @family sample functions
#' 
#' @inheritParams analyze_sampOptLHC
#' @inheritParams sample_srs
#'
#' @importFrom magrittr %>%
#' @importFrom methods is
#' 
#' @return An sf object with \code{n} stratified samples.
#' 
#' @export

sample_clhs <- function(mraster = mraster,
                        n = NULL,
                        existing = NULL,
                        access = NULL,
                        buff_inner = NULL,
                        buff_outer = NULL,
                        plot = FALSE,
                        details = FALSE,
                        ...)
{
  
  #--- Error management ---#
  
  if (!inherits(mraster, "SpatRaster"))
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  
  if (!is.numeric(n))
    stop("'n' must be type numeric")
  
  if (is.na(crs(mraster)))
    stop("'mraster' does not have a coordinate system")
  
  if (!is.logical(plot))
    stop("'plot' must be type logical")
  
  if (!is.logical(details))
    stop("'details' must be type logical")
  
  
  #--- determine crs of input mraster ---#
  crs <- crs(mraster)
  
  #--- access buffering if specified ---#
  
  if (!is.null(access)) {
    
    #--- error handling in the presence of 'access' ---#
    if (!inherits(access,"sf"))
      stop("'access' must be an 'sf' object")
    
    if(!inherits(sf::st_geometry(access),"sfc_MULTILINESTRING"))
      stop("'access' geometry type must be 'sfc_MULTILINESTRING'")
    
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
    access <- terra::vect(access)
    
    #--- make access buffer with user defined values ---#
    
    buff_in <- terra::buffer(x = access,
                             width = buff_inner)
    
    buff_out <- terra::buffer(x = access,
                              width = buff_outer)
    
    #--- make difference and aggregate inner and outer buffers to prevent sampling too close to access ---#
    
    buffer <- terra::aggregate(buff_out - buff_in)
    
    mraster_access <- terra::mask(mraster, mask = buffer)
    
    #--- extract covariate data from mraster ---#
    
    vals <- terra::as.data.frame(mraster_access, xy = TRUE, row.names = FALSE) %>%
      dplyr::rename(X = x,
                    Y = y)
    
  } else {
    
    #--- extract covariate data from mraster ---#
    
    vals <- terra::as.data.frame(mraster, xy = TRUE, row.names = FALSE) %>%
      dplyr::rename(X = x,
                    Y = y)
    
  }

  #--- Remove NA / NaN / Inf values ---#
  
  vals <- vals %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::mutate(type = "new")

  #--- error handing and existing samples preparation ---#
  
  if (!is.null(existing)){
    
    if (!inherits(existing, "data.frame") && !inherits(existing,"sf"))
      stop("'existing' must be a data.frame or sf object")
    
    #--- combined existing samples with vals dataframe ---#
    
    if(!inherits(existing,"sf")){
      
      if (any(! c("X", "Y") %in% colnames(existing)) ) {
        
        #--- if coordinate column names are lowercase change them to uppercase to match requirements ---#
        
        if (any(c("x", "y") %in% colnames(existing))) {
          
          existing <- existing %>%
            dplyr::rename(X = x,
                          Y = y)
          
          message("Column coordinates names for 'existing' are lowercase - converting to uppercase")
          
        } else {
          
          #--- if no x/y columns are present stop ---#  
          
          stop("'existing' must have columns named 'X' and 'Y'")
          
        }
        
      }
    
      existingSamples <- existing
    
    } else {
      
      existingSamples <- extract_metrics(mraster,existing,data.frame = TRUE)
      
    }
    
    #--- create dataset with labels for plotting ---#
    
    existingSamples <- existingSamples %>% 
        dplyr::mutate(type = "existing")
    
    #--- create conjoined existing dataset ---#
    
    vals <- rbind(existingSamples, vals)
    
  }
  
  vals_tp <- vals %>% dplyr::select(-type)
  
  ##########################
  #### SAMPLING ############
  ##########################
  
  #--- if existing samples are not provided ---#
  
  if (is.null(existing)){
    
    #--- remove 'type' during sampling ---#
    
    clhsOut <- clhs::clhs(x = vals_tp, size = n, ...)
      
    
      #--- if ... variables are provided the output is sometimes a list object ---#
    
      if(is.list(clhsOut)){
        
        samples <- clhsOut$sampled_data
        
      } else {
        
        #--- take samples from original vals dataframe for 'type' attribute ---#
        
        samples <- vals[clhsOut,]
        
      }
      
    } else {
      
      #--- same as above but this time including existing samples ---#

      clhsOut <- clhs::clhs(x = vals_tp, size = n, include = 1:nrow(existingSamples), ...)

      if( inherits(clhsOut, "list")){
        
        #--- locate row indices for each sample and extract ---#

        index <- clhsOut$index_samples
        
        samples <- vals[index,]

      } else {

        samples <- vals[clhsOut,]

      }

    }

  #--- convert coordinates to a spatial points object ---#
  samples <- samples %>%
    as.data.frame() %>%
    sf::st_as_sf(., coords = c("X", "Y"))
  
  #--- assign sraster crs to spatial points object ---#
  sf::st_crs(samples) <- crs
  
  if(isTRUE(plot)){
    
    if(!is.null(access)){
      
      #--- plot samples as well as aggregated access buffers ---#
      
      terra::plot(mraster[[1]])
      suppressWarnings(terra::plot(buffer, add = T, border = c("gray60"), col = "gray10", alpha = 0.1))
      suppressWarnings(terra::plot(samples, add = T, col = ifelse(samples$type=="existing","Black","Red")))
      
    } else {
    
    terra::plot(mraster[[1]])
    suppressWarnings(terra::plot(samples, add = T, col = ifelse(samples$type=="existing","Black","Red")))
    
    }
    
  }
  
  if ( isTRUE(details) ){
    
    #--- output metrics details along with stratification raster ---#
    
    output <- list(clhs = clhsOut, samples = samples)
    
    #--- output samples dataframe ---#
    return(output)
    
    
  } else {
    
    #--- just output raster ---#
    
    return(samples)
    
    
  }

  
}