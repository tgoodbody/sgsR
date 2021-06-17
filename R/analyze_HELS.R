#' Perform the Hypercube Evaluation of a Legacy Sample (HELS) algorithm using existing site data
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
#' 


analyze_HELS <- function(mraster = NULL,
                         existing = NULL,
                         nQuant = 20,
                         nSamp = 100,
                         retireSamp = TRUE,
                         plot = FALSE)
{
  
  if (!inherits(mraster,"SpatRaster"))
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  
  if (!inherits(existing, "data.frame") && !inherits(existing,"sf"))
    stop("'existing' must be a data.frame or sf object")
  
  if (!is.numeric(nQuant))
    stop("'nQuant' must be type numeric")
  
  if (!is.numeric(nSamp))
    stop("'nSamp' must be type numeric")
  
  if (!is.logical(retireSamp))
    stop("'retireSamp' must be type logical")
  
  #--- determine number of bands in 'mraster' ---#
  
  nb <- terra::nlyr(mraster)
  
  #--- determine crs of input sraster ---#
  crs <- crs(mraster)
  
  #--- extract covariates data from mraster ---#
  
  vals <- terra::as.data.frame(mraster, xy = TRUE, row.names = FALSE) %>%
    dplyr::rename(X = x,
           Y = y)
  
  #--- Remove NA / NaN / Inf values ---#
  
  vals <- vals %>%
    dplyr::filter(complete.cases(.))
  
  #--- Generate quantile matrix ---#
  
  mats <- analyze_popLHC(mraster = mraster,PCA = FALSE, nQuant = nQuant)
  
  #--- Change 0's to very small number to avoid division issues ---#
  
  mats$matCov[which(mats$matCov == 0)] <- 0.0000001
  
  #--- Create density matrix from covariates and length of mraster ---#
  
  matCovDens <- mats$matCov / nrow(vals)
  
  ###--- Prepare existing sample data ---###
  
  #--- extract covariates at existing sample locations ---#
  
  samples <- extract_metrics(mraster, existing, data.frame = TRUE)
  
  #--- Assign code to differentiate between original samples and those added during HELS algorithm ---#
  
  samples$type <- "existing"
  
  #--- Create data hypercube of existing samples to compare with mraster data ---#
  
  matCovSamp <- mat_cov(vals = samples[3:ncol(samples)], nQuant = nQuant, nb = nb, matQ = mats$matQ)
  
  #--- Change 0's to very small number to avoid division issues ---#
  
  matCovSamp[which(matCovSamp == 0)] <- 0.0000001
  
  #--- Create density matrix from covariates and length of mraster ---#
  
  matCovSampDens <- matCovSamp / nrow(samples)
  
  ###--- Selection of new samples based on density ---###
  
  #--- Ratio and ordering of data density and covariate density ---#
  
  ratio <- matCovSampDens / matCovDens
  
  #--- order the densities based on representation ---#
  
  #--- low to high ---#
  
  ratOrderUnder <- order(ratio)
  
  #--- high to low ---#
  
  ratOrderOver <- rev(ratOrderUnder)
  
  #--- Outline quantiles that are underrepresented (< 1) in the sample ---#
  
  underRep <- which(ratio < 1, arr.ind = TRUE)
  underRep <- cbind(underRep,which(ratio < 1))
  
  #--- Outline quantiles that are overrepresented (> 1) in the sample ---#
  
  overRep <- which(ratio > 1, arr.ind = TRUE)
  overRep <- cbind(overRep,which(ratio > 1))
  
  #--- begin sampling from highest discrepancy to lowest ---#
  
  newSamp <- nSamp
  position <- 1
  
  #--- begin while loop to sample ---#
  
  while(newSamp > 0){
  
    #--- determine the greatest discrepancy between sample and covariate data ---#
    
    repRankUnder <- which(underRep[,3] == ratOrderUnder[position])
    
    #--- determine row and column of most under represented quantile ---#
    
    repRow <- underRep[repRankUnder,1]
    repCol <- underRep[repRankUnder,2]
    
    #--- determine number of existing samples in selected quantile ---#
    
    sampExist <- floor(nrow(samples) * matCovSampDens[repRow,repcol])
    
    #--- determine max number of samples based on covariate density ---#
    
    sampOptim <- ceiling(nrow(samples) * matCovDens[repRow,repcol])
    
    #--- number of samples needed ---#
    
    sampNeed <- sampOptim - sampExist
  
    #--- we have a limited number of samples so we need to be sure not to over allocate ---#
    
    if(newSamp < sampNeed) sampNeed <- newSamp
    
    #--- selecting covariates based on quantile chosen ---#
    
    covLower <- mats$matQ[repRow,repCol]
    
    covUpper <- mats$matQ[repRow + 1,repCol]
    
    #--- subset covariate dataset for potential new samples ---#
    
    valsSub <- vals[vals[,(2 + repCol)] >= covLower & vals[,(2 + repCol)] <= covUpper,]
    
    #--- randomly sample within valsSub and extract randomly sampled cells ---#
    
    addSamp <- sample(nrow(valsSub), sampNeed)
    
    valsSubSamp <- valsSub[addSamp,]
    
    valsSubSamp$type <- "new"
    
    #--- remove samples from pool to ensure same cells are not sampled again ---#
    
    vals <- vals[-addSamp,]
    
    #--- add new samples to existing sample dataframe ---#
    
    samples <- rbind(samples,valsSubSamp)
    
    #--- update loop parameters ---#
    position <- position + 1
    newSamp <- newSamp - sampNeed
    print(sampNeed)

  }
  
  #--- Create data hypercube of existing samples to compare with mraster data ---#
  
  matCovSamp <- mat_cov(vals = samples[3:ncol(samples)], nQuant = nQuant, nb = nb, matQ = mats$matQ)
  
  #--- Change 0's to very small number to avoid division issues ---#
  
  matCovSamp[which(matCovSamp == 0)] <- 0.0000001
  
  #--- Create density matrix from covariates and length of mraster ---#
  
  matCovSampDens <- matCovSamp / nrow(samples)
  
  ###--- Selection of new samples based on density ---###
  
  #--- Ratio and ordering of data density and covariate density ---#
  
  ratio <- matCovSampDens / matCovDens
  
  #--- convert coordinates to a spatial points object ---#
  samples <- samples %>%
    as.data.frame() %>%
    sf::st_as_sf(., coords = c("X", "Y"))
  
  #--- assign sraster crs to spatial points object ---#
  sf::st_crs(samples) <- crs
  
  if(isTRUE(plot)){
    
    terra::plot(mraster[[1]])
    suppressWarnings(terra::plot(samples, add = T, col = ifelse(samples$type=="existing","Red","Black")))
    
  }
  
  return(ratio)
  
  
}