#' Evaluates existing sample data using an adapted Hypercube Evaluation of a Legacy Sample (HELS) algorithm
#' to remove redundant samples based on covariate and sample quantile comparisons
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
#' 

analyze_HELSrem <- function(mraster = NULL,
                            existing = NULL,
                            nQuant = 10
                            ){
  
  
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

#--- Remove quantiles that do not cover at least 1% area in eac covariate ---#

matCovDens[which(matCovDens <= 0.01)] <- NA

###--- Prepare existing sample data ---###

#--- extract covariates at existing sample locations ---#

samples <- extract_metrics(mraster, existing, data.frame = TRUE)

#--- Assign code to differentiate between original samples and those added during HELS algorithm ---#

samples$type <- "existing"
samples$n <- seq(1:nrow(samples))

#--- Rearrange columns ---#

samples <- samples %>%
  dplyr::select(X,Y,n,type,everything())

#--- Create data hypercube of existing samples to compare with mraster data ---#

matCovSamp <- mat_cov(vals = samples[5:ncol(samples)], nQuant = nQuant, nb = nb, matQ = mats$matQ)

#--- Change 0's to very small number to avoid division issues ---#

matCovSamp[which(matCovSamp == 0)] <- 0.0000001

#--- Create density matrix from covariates and length of mraster ---#

matCovSampDens <- matCovSamp / nrow(samples)

###--- Selection of new samples based on density ---###

#--- Ratio and ordering of data density and covariate density ---#

ratio1 <- matCovSampDens / matCovDens

#--- OUtline which quantiles are greater than 1 (overrepresented) ---#

ratioGr1 <- ratio1 > threshold

#--- set NA values to FALSE -- these are from quantiles with <= 1% of covariate coverage ---#

ratioGr1[is.na(ratioGr1)] <- FALSE

#--- Generate dataframe where quantiles that are TRUE in ratioGr1 are compared ---#
#--- indexOverRep will be used to interatively remove samples starting from the most overrepresented quantile in the sample data (Subject) ---#
#--- The subject will then be compared to overrepresented quantiles in all other covariates (Subject) to narrow down where samples should be prioritied to be removed ---#

#--- create empty dataframe to be populated ---#

indexOverRep <- data.frame()

#--- iterate through each row and column ---#

for(i in 1:nrow(ratioGr1)){
  for(j in 1:ncol(ratioGr1)){
    
    #--- SUBJECT - process when TRUE ---#
    
    if(isTRUE(ratioGr1[i,j])){
      
      #--- Do not compare quantiles from the same variable ---#
      
      for(ii in 1:nrow(ratioGr1)){
        for(jj in 1:ncol(ratioGr1)){
          
          #--- OBJECT - if the 2 quantiles share the same column - go next ---#
          
          if(jj == j) next
          
          #--- if the Subject == TRUE and Object == TRUE ---#
          
          if(isTRUE(ratioGr1[ii,jj])){
            
            #--- Determine overrepresented quantile range ---#
            
            covLower <- mats$matQ[ii, jj]
            covUpper <- mats$matQ[ii + 1, jj]
            
            #--- Populate dataframe listing the subject quantile (rowSub, colSub) and object quantile (rowObj, colObj) and associated OBJECT quantile ranges ---#
            
            df <- data.frame(rowSub = i, colSub = j, rowObj = ii, colObj = jj, covLower = covLower, covUpper = covUpper)
            
            #--- Rbind the dataframe together to create a comprensive comparative list of overrepresented quantiles ---#
            
            indexOverRep <- rbind(df,indexOverRep)
            
          } else {
            
            #--- if FALSE - go next ---#
            
            next
            
          }
          
        }
        
      }
      
    } else {
      
      #--- if FALSE - go next ---#
      
      next
      
    }
    
  }
  
}

#--- order the densities based on representation ---#
#--- low to high ---#

ratOrderUnder <- order(ratio1,na.last = NA)

#--- high to low ---#

ratOrderOver <- rev(ratOrderUnder)

#--- Outline quantiles that are overrepresented (> 1) in the sample ---#

overRep <- which(ratio1 > threshold, arr.ind = TRUE)
overRep <- cbind(overRep,which(ratio1 > threshold))

#--- determine decending order of overrepresentation ---#

indexOrder <- overRep %>%
  as.data.frame() %>%
  dplyr::rename(order = V3) %>%
  dplyr::arrange(match(order,ratOrderOver[1:nrow(overRep)])) %>%
  tidyr::unite(uniteObj, c("row","col"), remove = FALSE)

#--- looping parameters ---#

objRank <- 1
samplesNew <- samples

while(objRank <= nrow(overRep)){
  
  #--- determine row and column of most over represented quantile ---#
  
  repRow <- indexOrder[objRank,]$row
  repCol <- indexOrder[objRank,]$col
  
  #--- determine optimum number of samples based on covariate density ---#
  
  sampOptim <- ceiling(nrow(samples) * matCovDens[repRow,repCol])
  
  #--- select covariates based on quantile chosen ---#
  
  covLower <- mats$matQ[repRow,repCol]
  
  covUpper <- mats$matQ[repRow + 1,repCol]
  
  #--- create subject level subset of overrepresented samples ---#
  
  sub <- samplesNew %>% 
    dplyr::filter(.[,(4 + repCol)] >= covLower & .[,(4 + repCol)] <= covUpper)
  
  #--- Total number of sample matching subject quantile criteria ---#
  
  sampTot <- nrow(sub)
  
  #--- difference between optimal sample number and total samples meeting subject criteria ---#
  
  sampDiff <- abs(sampOptim - sampTot)
  
  #--- Subset indexOverRep to match subject criteria and arrange in order of descending overrepresentation of object quantiles ---#
  
  indexOverRepSub <- indexOverRep %>% 
    dplyr::filter(rowSub == repRow & colSub == repCol) %>%
    tidyr::unite(uniteObj, c("rowObj","colObj"), remove = FALSE) %>%
    dplyr::arrange(match(uniteObj,indexOrder$uniteObj))
  
  #--- loop interator to stop after subject / object pairings have been exhausted ---#
  
  subRank <- 1
  
  #--- Removal of samples loop ---#
  #--- Continue until ---#
  #--- 1) the sampDiff value has reached 0 (optimal sample number for subject quantile has been reached) ---#
  #--- 2) the subject / object pairings have been exhausted ---#
  #--- Samples that meet criteria will have $type set to "retired" ---#
  
  while(sampDiff > 0 && subRank <= nrow(indexOverRepSub)){
    
    for(i in 1:nrow(indexOverRepSub)){
      
      #--- update object order rank ---#
      
      subRank <- subRank + 1
      
      #--- create object level subset of overrepresented samples ---#
      
      subObj <- sub %>% 
        dplyr::filter(type != "retired") %>%
        dplyr::filter(.[,(4 + indexOverRepSub$colObj[i])] >= indexOverRepSub$covLower[i] & .[,(4 + indexOverRepSub$colObj[i])] <= indexOverRepSub$covUpper[i])
      
      #--- potential number of samples that can be retired ---#
      
      nRet <- nrow(subObj)
      
      #--- If there are no samples that match criteria - go next ---#
      
      if(nRet == 0) next
      
      #--- If the nRet > sampDiff - take the remaining number of samples to be removed ---#
      
      if(nRet > sampDiff){
        
        samplesNew <- subObj %>%
          dplyr::slice(sampDiff) %>%
          dplyr::mutate(type = "retired") %>%
          dplyr::rows_update(samplesNew, ., by= "n")
        
        #--- If sampDiff > nRet remove all samples available ---#
        
      } else {
        
        samplesNew <- subObj %>%
          dplyr::mutate(type = "retired") %>%
          dplyr::rows_update(samplesNew, ., by= "n")
        
      }
      
      #--- Update sampDiff iterator to reflect retirement of samples ---#
      
      sampDiff <- sampDiff - nRet
      
    }
    
  }
  
  #--- Update Subject order rank ---#
  
  objRank <- objRank + 1
  
}

#--- determine the difference im representation for quantiles ---#

#--- filter our all retired samples ---#
samples_ <- samplesNew %>% dplyr::filter(type != "retired")

matCovSamp_ <- mat_cov(vals = samples_[5:ncol(samples_)], nQuant = nQuant, nb = nb, matQ = mats$matQ)

#--- Change 0's to very small number to avoid division issues ---#

matCovSamp_[which(matCovSamp_ == 0)] <- 0.0000001

#--- Create density matrix from covariates and length of mraster ---#

matCovSampDens_ <- matCovSamp_ / nrow(samples_)

###--- Selection of new samples based on density ---###

#--- Ratio and ordering of data density and covariate density ---#

ratio1_ <- matCovSampDens_ / matCovDens

rdiff <- ratio1_ - ratio1

# par(mfrow=c(1,3))
# par(mar=c(5.1, 4.1, 4.1, 4.1))
# 
# plot(ratio1, breaks=range(ratio1_,na.rm=TRUE), col=topo.colors, main = "Original Data", asp=TRUE, axis.col=NULL, axis.row=NULL, xlab='', ylab='')
# plot(ratio1_, breaks=range(ratio1_,na.rm=TRUE), col=topo.colors, main = "Reduced Data", asp=TRUE, axis.col=NULL, axis.row=NULL, xlab='', ylab='')
# plot(rdiff, breaks=range(rdiff,na.rm=TRUE), col=topo.colors, main = "Difference", asp=TRUE, axis.col=NULL, axis.row=NULL, xlab='', ylab='')

k <- reshape2::melt(ratio1, c("y", "x"), value.name = "z") %>%
  dplyr::filter(complete.cases(.))

k$type <- "pre"

k1 <- reshape2::melt(ratio1_, c("y", "x"), value.name = "z") %>%
  dplyr::filter(complete.cases(.)) 

k1$type <- "post"

k2 <- reshape2::melt(rdiff, c("y", "x"), value.name = "z") %>%
  dplyr::filter(complete.cases(.))

k2$type <- "diff"


kk <- rbind(k,k1)

ggplot(kk, aes(x=x, y=y)) + geom_tile(aes(fill = z), colour = "grey50") +
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0.7) + facet_grid(.~type)


p1 <- ggplot(k, aes(x=x, y=y, fill = z)) + geom_tile(colour = "grey50") + scale_fill_viridis() + theme_light()
p2 <- ggplot(k1, aes(x=x, y=y, fill = z)) + geom_tile(colour = "grey50") + scale_fill_viridis() + theme_light()
p3 <- ggplot(k2, aes(x=x, y=y, fill = z)) + geom_tile(colour = "grey50") + scale_fill_viridis() + theme_light()


ggplot(k, aes(x=x, y=y, fill = z)) + geom_tile(colour = "grey50") + scale_colour_gradientn(colours = c("red","white","darkblue"), values = c(min(k$z), median(k$z), max(k$z)))


gridExtra::grid.arrange(p1, p2, p3, nrow = 1, ncol = 3)

#--- convert coordinates to a spatial points object ---#
samples <- samplesNew %>%
  as.data.frame() %>%
  sf::st_as_sf(., coords = c("X", "Y"))

#--- assign sraster crs to spatial points object ---#
sf::st_crs(samples) <- crs



terra::plot(mraster[[1]])
terra::plot(samples, add = T, col = ifelse(samples$type=="existing","Red","Black"))






########################################





print(cutTotal)

matCovSamp <- mat_cov(vals = samples[5:ncol(samples)], nQuant = nQuant, nb = nb, matQ = mats$matQ)

#--- Change 0's to very small number to avoid division issues ---#

matCovSamp[which(matCovSamp == 0)] <- 0.0000001

#--- Create density matrix from covariates and length of mraster ---#

matCovSampDens <- matCovSamp / nrow(samples)

###--- Selection of new samples based on density ---###

#--- Ratio and ordering of data density and covariate density ---#

ratio2 <- matCovSampDens / matCovDens
print(ratio1 - ratio2)




# model <- prcomp(vals[,3:5], scale. = TRUE, center = TRUE)
# w=TRUE)
# scores_all <- data.frame(model$x[,1:2])
# scores_all$type <- "all"
# 
# scores_samp <- data.frame(predict(model,samples[,5:7]))[,1:2]
# scores_samp$type <- "sample"
# 
# df <- rbind(scores_all,scores_samp)
# 
# ggplot(df,aes(x=PC1,y=PC2,colour = type)) + geom_point(size =2,) + labs(title="Plotting Customer Data against PC1 and PC2")
# 
# 
# 
# plot(scores$PC1,scores$PC2)
# par(ne
#     plot(new$PC1,new$PC2,col="red",add=T)

}
    
    