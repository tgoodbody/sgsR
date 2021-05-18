#' Population level analysis of metric raster data to determine optimal Latin Hypercube sample size
#' @family analyze functions
#' 
#' @inheritParams analyze_popLHC
#' 
#' @param minSamp Numeric. Minimum sample size to test. \code{default = 10}.
#' @param maxSamp Numeric. Maximum sample size to test. \code{default = 100}.
#' @param step Numeric. Sample step size for each iteration. \code{default = 10}.
#' @param rep Numeric. Internal repetitions for each sample size. \code{default = 10}.
#' @param iter Numeric. Internal to \code{clhs} - A positive number, giving the number of 
#' iterations for the Metropolis-Hastingsannealing process. Defaults to \code{10000}.
#' 
#' @importFrom methods is
#' 
#' @return TBD
#' 
#' @export


analyze_sampOptLHC <- function(popLHC = NULL,
                               PCA = TRUE,
                               quant = TRUE,
                               KLdiv = TRUE,
                               minSamp = 10,
                               maxSamp = 100,
                               step = 10,
                               rep = 10,
                               iter = 10000)
{
  
  if (!is.list(popLHC))
    stop("'popLHC' must be a list")
  
  if (any(!names(popLHC) %in% c("values", "pcaLoad", "matQ","matCov" )))
    stop(paste0("'popLHC' must be the output from the 'analyze_popLHC()' function"))
  
  if (!is.logical(PCA))
    stop("'PCA' must be type logical")
  
  if (!is.logical(quant))
    stop("'quantiles' must be type logical")
  
  if (!is.logical(KLdiv))
    stop("'KLdiv' must be type logical")
  
  if (!is.numeric(minSamp))
    stop("'minSamp' must be type numeric")
  
  if (!is.numeric(maxSamp))
    stop("'maxSamp' must be type numeric")
  
  if (!is.numeric(step))
    stop("'step' must be type numeric")
  
  
  #--- set global variables ---#
  
  nb <- ncol(popLHC$values)
  
  #--- Establish sampling sequence ---#
  
  sampSeq <- seq(minSamp,maxSamp,step)
  
  matSeq <- matrix(NA, ncol = 6, nrow = length(sampSeq))
  
  #--- Apply functions for each potential sample size ---#
  
  for (tSamp in 1:length(sampSeq)){ 
    
    P <- matrix(NA, ncol=6, nrow=length(sampSeq)) # placement for iteration outputs
    
    #--- run conditionel latin hypercube sampling ---#
    
    for (j in 1:rep){ #Note that this takes quite a while to run to completion
        
        #--- perform conditionel Latin Hypercube Sampling ---#
        
        ss <- clhs::clhs(popLHC$values, size = sampSeq[tSamp], progress = TRUE, iter = iter) 
        
        samples <- popLHC$values[ss,]
      
      
      # --- PCA similarity factor testing ---#
      
      if(isTRUE(PCA)){
        
        #--- perform PCA analysis for the samples to determine variance in each component ---#
        
        pcaS <- stats::prcomp(samples, scale=TRUE, center=TRUE)
        
        #--- extract sample pca scores ---#
        
        pcaSScores <-  as.data.frame(pcaS$x)
        
        #--- extract sample PCA loadings ---#
        
        pcaLoadSamp <- matrix(NA, ncol=nb, nrow=nb)
        
        for (i in 1:nb){
          
          pcaLoadSamp[i,]<- as.matrix(t(pcaS$rotation[i,]))
          
        }
        
        #--- Perfrom the Krznowski 1979 calculation ---#
        
        pop <- popLHC$pcaLoad[,1:2]  
        samp <- pcaLoadSamp[,1:2]
        
        #--- transpose matrices ---#
        
        popT <- t(pop)
        sampT <- t(samp)
        
        #--- matrix multiplication ---#
        
        S <- popT %*% samp %*% sampT %*% pop
        
        matFinal[j,1] <- sum(diag(S)) / 2
        
        
      }
      
      #--- Comparison of quantiles ---#
      
      if(isTRUE(quant)){
        
        for(var in 1:nb){
          
          #--- Calculate sample quantiles ---#
          
          sampleQuant <- quantile(samples[,var], probs = seq(0, 1, 0.25), names = F, type = 7)
          
          #--- Calculate population quantiles ---#
          
          popQuant <- quantile(popLHC$values[,var], probs = seq(0, 1, 0.25),names = F, type = 7)
          
          #--- populate quantile differences into matFinal ---#
          
          matFinal[j,var+1] <- sqrt((popQuant[1]-sampleQuant[1])^2 +
                                      (popQuant[2]-sampleQuant[2])^2 +
                                      (popQuant[3]-sampleQuant[3])^2 +
                                      (popQuant[4]-sampleQuant[4])^2 )
        }
        
        #--- calculate mean distance from all variables calculated above ---#
        
        matFinal[j,6] <- mean(matFinal[j,2:sum((1+nb))])
        
      }
      
      if(isTRUE(KLdiv)){
        
        #--- calculate sample covariate matrix ---#
        
        sampleCov <- mat_cov(vals = samples,
                             nQuant = nrow(popLHC$matCov),
                             nb = nb,
                             matQ = popLHC$matQ)
        
        #--- calculate KL divergence ---#
        
        #--- create empty vector for populating in loop ---#
        
        kld.vars <- c()
        
        #--- Generate KL divergence for each metric and populate kld.vars ---#
        
        for(kl in 1:nb){
          
          #--- calculate divergence ---#  
          
          kld <- entropy::KL.empirical(popLHC$matCov[,kl], sampleCov[,kl])
          
          #--- populate vector ---#
          
          kld.vars[kl] <- kld
          
        }
        
        #--- calculate mean divergence ---#
        #--- Divergence of 0 means samples do not diverge from pop ---#
        
        KLout <- mean(kld.vars)
        
        #--- populate to final matrix ---#
        
        matFinal[j,7]<- KLout  
        
      }
      
    }
    
    #--- create outputs for all tests ---#
    
    #arrange outputs
    matSeq[tSamp,1] <- mean(matFinal[,6])
    matSeq[tSamp,2] <- sd(matFinal[,6])
    matSeq[tSamp,3] <- min(matFinal[,1])
    matSeq[tSamp,4] <- max(matFinal[,1])
    matSeq[tSamp,5] <- mean(matFinal[,7])
    matSeq[tSamp,6] <- sd(matFinal[,7])
    
    
  }
  
  dat.seq<- as.data.frame(cbind(sampSeq,matSeq))
  names(dat.seq)<- c("samp_nos", "mean_dist","sd_dist", "min_S", "max_S", "mean_KL","sd_KL")
  
  
  return(dat.seq)
}

