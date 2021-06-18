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
  select(X,Y,n,type,everything())

#--- Create data hypercube of existing samples to compare with mraster data ---#

matCovSamp <- mat_cov(vals = samples[5:ncol(samples)], nQuant = nQuant, nb = nb, matQ = mats$matQ)

#--- Change 0's to very small number to avoid division issues ---#

matCovSamp[which(matCovSamp == 0)] <- 0.0000001

#--- Create density matrix from covariates and length of mraster ---#

matCovSampDens <- matCovSamp / nrow(samples)

###--- Selection of new samples based on density ---###

#--- Ratio and ordering of data density and covariate density ---#

ratio1 <- matCovSampDens / matCovDens

#--- We want to make sure we dont take away
# kk <- rowSums((ratio1 < 1)/ncol(ratio1), na.rm = TRUE)
# 
# kk[which(kk == 0)] <- 0.01
# 
# k <- ratio1 * kk

print(ratio1)

ratioGr1 <- ratio1 > 1

ratioGr1[is.na(ratioGr1)] <- FALSE

test <- data.frame()
for(i in 1:nrow(ratioGr1)){
  for(j in 1:ncol(ratioGr1)){
    
    if(isTRUE(ratioGr1[i,j])){
      
      for(ii in 1:nrow(ratioGr1)){
        for(jj in 1:ncol(ratioGr1)){
          
          if(jj == j) next
          
          if(isTRUE(ratioGr1[ii,jj])){
            
            covLower <- mats$matQ[ii, jj]
            covUpper <- mats$matQ[ii + 1, jj]
            
            df <- data.frame(rowSub = i, colSub = j, rowObj = ii, colObj = jj, covLower = covLower, covUpper = covUpper)
            
            test <- rbind(df,test)
          } else {
            
            next
            
          }
          
          
        }
      }
      
      
    } else {
      
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

overRep <- which(ratio1 > 1, arr.ind = TRUE)
overRep <- cbind(overRep,which(ratio1 > 1))

xx <- overRep %>%
  as.data.frame() %>%
  arrange(match(V3,ratOrderOver[1:nrow(overRep)])) %>%
  unite(uniteObj, c("row","col"), remove = FALSE)

cutTotal <- 0
position <- 1

while(position <= nrow(overRep)){
  
  #--- determine row and column of most under represented quantile ---#
  
  repRow <- xx[position,]$row
  repCol <- xx[position,]$col
  
  #--- determine max number of samples based on covariate density ---#
  
  sampOptim <- ceiling(nrow(samples) * matCovDens[repRow,repcol])
  
  #--- selecting covariates based on quantile chosen ---#
  
  covLower <- mats$matQ[repRow,repCol]
  
  covUpper <- mats$matQ[repRow + 1,repCol]
  
  #--- subset covariate dataset for potential new samples ---#
  set.seed(420)
  sub <- samples %>% 
    dplyr::filter(samples[,(4 + repCol)] >= covLower & samples[,(4 + repCol)] <= covUpper)
  
  t1 <- test %>% 
    filter(rowSub == repRow & colSub == repCol) %>%
    unite(uniteObj, c("rowObj","colObj"), remove = FALSE) %>%
    arrange(match(uniteObj,xx$uniteObj))
  
  
  while(length(sub) > length(sampOptim)){
    
    for(i in 1:nrow(t1)){
      
      subsub <- sub %>% 
        filter(sub[,(4 + t1$colObj[i])] >= t1$covLower[i] & sub[,(4 + t1$colObj[i])] <= t1$covUpper[i])
      
      if(nrow(subsub) == 0) next
      
      
    }
  }
  
  sub 
  
  cut <- abs(sampOptim-nrow(sub))
  print(cut)
  
  samples <- sub %>%
    dplyr::sample_n(cut) %>%
    dplyr::mutate(type = "retired") %>%
    dplyr::rows_update(samples, ., by= "n")
  
  samples <- samples %>% filter(type != "retired")
  
  cutTotal <- cutTotal + cut
  
  position <- position + 1
  
}

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


