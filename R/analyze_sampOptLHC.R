#' Analyze optimal LHC sample number
#' 
#' @description Population level analysis of metric raster data to determine optimal Latin Hypercube sample size
#' @family analyze functions
#' 
#' @inheritParams analyze_popLHC
#' 
#' @param popLHC List. Output from \code{analyze_popLHC} function.
#' @param minSamp Numeric. Minimum sample size to test. \code{default = 10}.
#' @param maxSamp Numeric. Maximum sample size to test. \code{default = 100}.
#' @param step Numeric. Sample step size for each iteration. \code{default = 10}.
#' @param rep Numeric. Internal repetitions for each sample size. \code{default = 10}.
#' @param iter Numeric. Internal to \code{clhs} - A positive number, giving the number of 
#' iterations for the Metropolis-Hastings annealing process. Defaults to \code{10000}.
#' 
#' @importFrom methods is
#' 
#' @return data.frame with summary statistics.
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
  
  #--- Set global vars ---#
  x <- y <- NULL
  
  #--- Error handling ---#
  
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
  
  matSeq <- matrix(NA, nrow = length(sampSeq), ncol = 6)
  
  #--- Apply functions for each potential sample size ---#
  
  for (tSamp in 1:length(sampSeq)){ 
    
    #--- matrix holding iteration outputs ---#
    
    matFinal <- matrix(NA, nrow=rep, ncol=7) 
    
    #--- run conditional latin hypercube sampling ---#
    
    for (j in 1:rep){
        
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
          
          sampleQuant <- stats::quantile(samples[,var], probs = seq(0, 1, 0.25), names = F, type = 7)
          
          #--- Calculate population quantiles ---#
          
          popQuant <- stats::quantile(popLHC$values[,var], probs = seq(0, 1, 0.25),names = F, type = 7)
          
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

    matSeq[tSamp,1] <- mean(matFinal[,6])
    matSeq[tSamp,2] <- stats::sd(matFinal[,6])
    matSeq[tSamp,3] <- min(matFinal[,1])
    matSeq[tSamp,4] <- max(matFinal[,1])
    matSeq[tSamp,5] <- mean(matFinal[,7])
    matSeq[tSamp,6] <- stats::sd(matFinal[,7])
    
    
  }
  
  dfFinal<- as.data.frame(cbind(sampSeq,matSeq))
  names(dfFinal)<- c("n", "mean_dist","sd_dist", "min_S", "max_S", "mean_KL","sd_KL")
  
  #--- plot the outputs and determine optimal sample size based on mean_KL divergence ---#
  
  plot_LHCOptim(dfFinal,
                maxSamp)
  
  return(dfFinal)
}


plot_LHCOptim <- function(dfFinal,
                          maxSamp){
  
  #--- set global vars ---#
  
  df.x <- NULL
  
  #--- Create normalized ---#
  df <- data.frame(x = dfFinal[,1],
                        y = 1-(dfFinal[,6]-min(dfFinal[,6]))/(max(dfFinal[,6])-min(dfFinal[,6])))
  
  #Parametise Exponential decay function
  plot(df$x, df$y, xlab="sample number", ylab= "1 - KL Divergence")          # Initial plot of the data
  
  # Prepare a good inital state
  theta.0 <- max(df$y) * 1.1
  model.0 <- stats::lm(log(- y + theta.0) ~ x, data=df)
  alpha.0 <- -exp(stats::coef(model.0)[1])
  beta.0 <- stats::coef(model.0)[2]
  
  start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
  
  # Fit the model
  model <- stats::nls(y ~ alpha * exp(beta * x) + theta , data = df, start = start)
  
  
  # add fitted curve
  
  predicted <- stats::predict(model, list(x = df$x))
  
  plot(df$x, df$y, xlab = "# of samples", ylab = "norm mean KL divergence")
  graphics::lines(df$x, predicted, col = 'skyblue', lwd = 3)
  
  x1<- c(-1, maxSamp)
  y1<- c(0.95, 0.95)
  graphics::lines(x1,y1, lwd=2, col="red")
  
  num <- data.frame(df$x,predicted) %>% 
    dplyr::filter(abs(predicted - 0.95) == min(abs(predicted - 0.95))) %>% 
    dplyr::select(df.x) %>%
    dplyr::pull()
  
  message(paste0("Your optimum estimated sample size based on KL divergence is: ",num))
  
  x2<- c(num, num); y2<- c(0, 1)
  graphics::lines(x2,y2, lwd=2, col="red")
  
}

