#' Perform the COunt of ObServations (COOBS) algorithm using existing site data
#' and raster metrics
#' @family analyze functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams extract_existing
#' 
#' @param threshold Numeric. Proxy maximum pixel quantile to avoid outliers. \code{default = 0.95}
#' 
#' @return 
#' 
#' @importFrom magrittr %>%
#' @importFrom methods is
#' @importFrom foreach %dopar%
#' 
#' @export


analyze_COOBS <- function(mraster = NULL,
                          existing = NULL,
                          cores = 1,
                          threshold = 0.95,
                          plot = FALSE)
  {
  
  #--- Error handling ---#

  if (!inherits(mraster,"SpatRaster"))
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  
  if (!inherits(existing, "data.frame") && !inherits(existing,"sf"))
    stop("'existing' must be a data.frame or sf object")
  
  if (!is.numeric(threshold))
    stop("'threshold' must be type numeric")
  
  #--- extract covariates data from mraster ---#
  
  vals <- terra::as.data.frame(mraster, xy = TRUE, row.names = FALSE)
  
  #--- Remove NA / NaN / Inf values ---#
  
  vals <- vals %>%
    dplyr::filter(complete.cases(.))
  
  #--- Generate covariance matrix ---#
  
  covMat <- as.matrix(cov(vals[,3:ncol(vals)]))
  
  #--- extract covariates at existing sample locations ---#
  
  samples <- extract_metrics(mraster, existing, data.frame = TRUE)
  
  #--- create parallel processing structure ---#
  
  cl <- snow::makeCluster(spec = cores)
  doParallel::registerDoParallel(cl)
  
  #--- iterate parallel processing of mahalanobis distance ---#
  
  loop <- foreach::foreach(i = 1:nrow(vals), .combine = "c") %dopar% {
    
    cell <- vals[i,3:ncol(vals)]
    
    #--- Determine distance for each pixel in raster ---#
    
    pixDist <- stats::mahalanobis(x = as.matrix(vals[,3:ncol(vals)]), center = as.matrix(cell), cov = covMat)
    
    #--- Determine min and max distance values for each ---#
    
    pixMin <- min(pixDist)
    pixMax <- quantile(pixDist, probs = threshold)
    
    #--- Determine distance for each sample location ---#
    
    sampDist<- mahalanobis(x = as.matrix(samples[,3:ncol(samples)]), center = as.matrix(cell), cov = covMat) #calculate distance of observations to all other pixels
    
    #--- Normalize distance between data and samples)
    
    sampNDist <- (sampDist - pixMin) / (pixMax - pixMin) 
    
    #--- If sampDist > 1 sampDist > maxDist ---#
    
    sampNDist[sampNDist > 1] <- 1 
    
    #--- larger values equate to more similarity ---#
    
    sampNDist <- 1 - sampNDist
    
    #--- establish count above threshold ---#
    
    sum(sampNDist >= threshold)

    
  }
  
  #--- End parallel ---# 
  snow::stopCluster(cl)
  
  #--- Coerce output from parallel to a new attribute in covariates ---#
   
  vals$nSamp <- loop
  
  #--- convert nSamp to raster ---#
  
  rout <- terra::rast(as.matrix(vals[,c("x", "y", "nSamp")]), type = "xyz")
  
  return(rout)
  
  #--- Plot output ---#
  
  if (isTRUE(plot)){
  
    plot(r1)
    plot(existing, add = TRUE)
    
  }

}