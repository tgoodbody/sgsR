#' Analyze covariates for LHC
#' 
#' @description Population level analysis of metric raster data to determine optimal Latin Hypercube sample size
#' @family analyze functions
#'
#' @inheritParams strat_kmeans
#' 
#' @param PCA Logical. Calculates principal component loadings of the population for PCA similarity factor testing. 
#' \code{default = TRUE}.
#' @param quant Logical. Calculates quantile matrix of the population for quantile comparison testing. 
#' \code{default = TRUE}.
#' @param nQuant Numeric. Number of quantiles to divide the population into for \code{quant}. 
#' \code{default = 20}.
#' @param KLdiv Logical. Calculates covariate matrix of the population for Kullback–Leibler divergence testing. 
#' \code{default = TRUE}. Relies on \code{quant = TRUE} to calculate.
#' 
#' @importFrom methods is
#' 
#' 
#' @return List of matrices to be used as input for \code{analyze_sampOptLHC}.
#' 
#' @export

analyze_popLHC <- function(mraster,
                        PCA = TRUE,
                        quant = TRUE,
                        nQuant = 20,
                        KLdiv = TRUE)
{
  
  if (!inherits(mraster,"SpatRaster"))
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  
  if (!is.logical(PCA))
    stop("'PCA' must be type logical")
  
  if (!is.logical(quant))
    stop("'quantiles' must be type logical")
  
  if (!is.logical(KLdiv))
    stop("'KLdiv' must be type logical")
  
  #--- determine number of bands in 'mraster' ---#
  
  nb <- terra::nlyr(mraster)
  
  #--- Extract values from mraster ---#
  
  vals <- terra::values(mraster)
  
  #--- Determine index of each cell so to map values correctly without NA ---#
  
  vals[!is.finite(vals)] <- NA
  
  #--- Remove NA / NaN / Inf values ---#
  
  vals <- vals %>%
    as.data.frame() %>%
    dplyr::filter(stats::complete.cases(.))
  
  #--- PCA loadings for the population ---#
  
  if(isTRUE(PCA)){
    
    #--- perform PCA analysis for the population to determine variance in each component ---#
    
    pca <- stats::prcomp(vals, scale=TRUE, center=TRUE)
    
    #--- extract pca scores ---#
    
    pcaScores <-  as.data.frame(pca$x)
    
    #--- extract PCA loadings ---#
    
    pcaLoad <- matrix(NA, ncol=nb, nrow=nb)
    
    for (i in 1:nb){
      
      pcaLoad[i,] <- as.matrix(t(pca$rotation[i,]))
      
    }
    
  } else {
    
    pcaLoad <- NULL
    
  }
  
  #--- Quantiles of the population ---#
  
  if(isTRUE(quant)){
    
    if (!is.numeric(nQuant))
      stop("'nQuantiles' must be type numeric")
    
    matQ <- mat_quant(vals,
                         nQuant,
                         nb)
    
  } else {
    
    matQ <- NULL
    
  }
  
  #--- covariate hypercube for KL divergence test ---#
  
  if(isTRUE(KLdiv)){
    
    if(is.null(matQ))
      stop("KL divergence requires quantile matrix creation. Set 'quant == TRUE'.")
    
    #--- create covariate matrix of the quantiles ---#
    
    message("Creating covariance matrix.")
    
    matCov <- mat_cov(vals, nQuant, nb, matQ)
    
  } else {
    
    matCov <- NULL
    
  }
  
  #--- create list of outputs ---#
  
  lout <- list(values = vals,
               pcaLoad = pcaLoad,
               matQ = matQ,
               matCov = matCov)
  
  #--- remove NULL list objects ---#
  
  lout <- lout[!sapply(lout,is.null)]
  
  #--- output list ---#
  
  return(lout)
                       
}

#--- function for calculating quantile matrix ---#

mat_quant <- function(vals,
                      nQuant,
                      nb)
{
  
  #--- create covariance matrix of the quantiles ---#
  
  matQ <- matrix(NA, nrow =(nQuant+1), ncol = nb)
  
  for(i in 1:nb){
    
    #--- calculate min, max, and range of vals for covariates ---#
    
    maxVal <- max(vals[,i])
    minVal <- min(vals[,i])
    rangeVal <- maxVal - minVal
    
    #--- determine the width of each quantile ---#
    
    widthQuant <- rangeVal/nQuant 
    
    #--- populate quantile matrix ---#
    
    matQ[,i]<- seq(from = minVal, to = maxVal, by = widthQuant)
    
  }
  
  return(matQ)
  
}

#--- function for calculating covariance matrix with a progress bar---#

mat_cov <- function(vals,
                    nQuant,
                    nb,
                    matQ)
{
  
  matCov <- matrix(0, nrow = nQuant, ncol = nb)
  
  #--- create progress text bar ---#
  
  iterations <- nrow(vals)
  pb <- utils::txtProgressBar(min = 1, max = iterations, style = 3)
  
  #--- for each row in dataframe ---#
  
  for (i in 1:nrow(vals)){
    
    setTxtProgressBar(pb,i)
    
    count <- 1 
    
    #--- for each column in dataframe ---#
    
    for (j in 1:nb){ 
      
      indv <- vals[i,j] 
      
      #--- for each quantile in covariates ---#
      
      for (k in 1:nQuant){  
        
        #--- determine upper and lower limits ---#
        
        limL <- matQ[k, count] 
        
        limU <- matQ[k+1, count] 
        
        #--- determine whether the sample row fits within the quantile being analyzed ---#
        
        if (indv >= limL & indv <= limU){
          
          matCov[k, count] <- matCov[k, count] + 1
          
        }
        
      }
      
      count <- count + 1
      
    }
    
  }
  
  close(pb)
  
  return(matCov)
  
}  

#--- function for calculating covariance matrix without a progress bar---#

mat_covNB <- function(vals,
                    nQuant,
                    nb,
                    matQ)
{
  
  matCov <- matrix(0, nrow = nQuant, ncol = nb)
  
  #--- for each row in dataframe ---#
  
  for (i in 1:nrow(vals)){

    
    count <- 1 
    
    #--- for each column in dataframe ---#
    
    for (j in 1:nb){ 
      
      indv <- vals[i,j] 
      
      #--- for each quantile in covariates ---#
      
      for (k in 1:nQuant){  
        
        #--- determine upper and lower limits ---#
        
        limL <- matQ[k, count] 
        
        limU <- matQ[k+1, count] 
        
        #--- determine whether the sample row fits within the quantile being analyzed ---#
        
        if (indv >= limL & indv <= limU){
          
          matCov[k, count] <- matCov[k, count] + 1
          
        }
        
      }
      
      count <- count + 1
      
    }
    
  }
  
  return(matCov)
  
}


