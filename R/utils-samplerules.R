#' Sampling rules
#'
#' @inheritParams sample_strat
#' @inheritParams sample_ahels
#' @family rules
#' @name rules
NULL

#' Stratified rule 1
#' @family rules
#' @rdname rules
#' @keywords internal
#' @export

strat_rule1 <- function(n, #number of samples
                        s, #strata number being sampled
                        i, #loop iteration
                        strat_mask,
                        add_strata,
                        extraCols,
                        mindist){
  
  #--- create indices for all, NA, and valid sampling candidates ---#
  
  idx_all <- 1:terra::ncell(strat_mask)
  idx_na <- !complete.cases(terra::values(strat_mask))
  validCandidates <- idx_all[!idx_na]
  
  #--- Rule 1 sampling ---#
  nCount <- 0 # Number of sampled cells
  
  #--- if areas exist where cells are grouped to meet moving window standard ---#
  
  if(length(validCandidates != 0)){
    
    # While loop for RULE 1
    while (length(validCandidates) > 0 & nCount < n) {
      #-- Identify potential sample from candidates ---#
      smp <- sample(1:length(validCandidates), size = 1)
      
      smp_cell <- validCandidates[smp]
      
      #--- Remove sampled cell from validCandidates so that it cannot be sampled again later ---#
      validCandidates <- validCandidates[-smp]
      
      #--- extract coordinates and sample details ---#
      
      add_temp <- data.frame(
        cell = smp_cell,
        X = terra::xFromCell(strat_mask, smp_cell),
        Y = terra::yFromCell(strat_mask, smp_cell),
        strata = s
      )
      
      #--- populate add_temp with values ---#
      add_temp$type <- "new"
      add_temp$rule <- "rule1"
      add_temp[, extraCols] <- NA
      
      #--- If add_strata is empty, sampled cell accepted ---#
      
      if (nrow(add_strata) == 0) {
        add_strata <- add_temp[, c("cell", "X", "Y", "strata", "type", "rule", extraCols)]
        
        nCount <- nCount + 1
        
      } else {
        
        #--- If add_strata isn't empty, check distance with all other sampled cells in strata ---#
        
        if (!is.null(mindist)) {
          dist <- spatstat.geom::crossdist(add_temp$X, add_temp$Y, add_strata$X, add_strata$Y)
          
          #--- If all less than 'mindist' - accept sampled cell otherwise reject ---#
          if (all(as.numeric(dist) > mindist)) {
            add_strata <- rbind(add_strata, add_temp[, c("cell", "X", "Y", "strata", "type", "rule", extraCols)])
            
            nCount <- nCount + 1
          }
        } else {
          
          #--- if mindist is not defined ---#
          
          if (add_temp$cell %in% add_strata$cell)  next
          
          add_strata <- rbind(add_strata, add_temp[, c("cell", "X", "Y", "strata", "type", "rule", extraCols)])
          
          nCount <- nCount + 1
        }
      }
      
      return(list(add_strata = add_strata, nCount = nCount))
    }
    
  } else if (i == 1) {
    
    add_strata <- data.frame(
      cell = NULL,
      X = NULL,
      Y = NULL,
      strata = NULL,
      type = NULL,
      rule = NULL
    )
    
    #--- populate add_temp with values ---#
    add_strata[, extraCols] <- NULL
    
    #--- if no areas exist where cells are grouped to meet moving window standard ---#
    
    return(list(add_strata = add_strata, nCount = nCount))
    
  } else if (i != 1) {
    
    add_strata <- data.frame(
      cell = NULL,
      X = NULL,
      Y = NULL,
      strata = NULL,
      type = NULL,
      rule = NULL
    )
    
    #--- populate add_temp with values ---#
    add_strata[, extraCols] <- NULL
    
    return(list(add_strata = add_strata, nCount = nCount))
    
  }
  
}

#' Stratified rule 2
#' @family rules
#' @rdname rules
#' @keywords internal
#' @export

strat_rule2 <- function(n, #number of samples
                        s, #strata number being samples
                        add_strata, #sampled output form rule 1
                        nCount, #number of allocated samples from rule 1
                        strata_m,
                        extraCols,
                        mindist){
  
  ### --- RULE 2 sampling ---###
  
  if (nCount < n) {
    idx_all <- 1:terra::ncell(strata_m)
    idx_na <- !complete.cases(terra::values(strata_m))
    validCandidates <- idx_all[!idx_na]
    
    while (length(validCandidates) > 0 & nCount < n) {
      
      #-- identify potential sample from candidates ---#
      smp <- sample(1:length(validCandidates), size = 1)
      
      smp_cell <- validCandidates[smp]
      
      #--- Remove sampled cell from validCandidates so that it cannot be sampled again later ---#
      
      validCandidates <- validCandidates[-smp]
      
      #--- extract coordinates and sample details ---#
      
      add_temp <- data.frame(
        cell = smp_cell,
        X = terra::xFromCell(strata_m, smp_cell),
        Y = terra::yFromCell(strata_m, smp_cell),
        strata = validCandidates[smp_cell]
      )
      
      add_temp$rule <- "rule2"
      add_temp$type <- "new"
      add_temp[, extraCols] <- NA
      add_temp$strata <- s
      
      if (nrow(add_strata) == 0) {
        add_strata <- add_temp[, c("cell", "X", "Y", "strata", "type", "rule", extraCols)]
        
        nCount <- nCount + 1
        
      } else {
      
        #--- If add_strata isn't empty, check distance with all other sampled cells in strata ---#
        
        if (!is.null(mindist)) {
          dist <- spatstat.geom::crossdist(add_temp$X, add_temp$Y, add_strata$X, add_strata$Y)
          
          #--- If all less than 'mindist' - accept sampled cell otherwise reject ---#
          if (all(as.numeric(dist) > mindist)) {
            add_strata <- rbind(add_strata, add_temp[, c("cell", "X", "Y", "strata", "type", "rule", extraCols)])
            
            nCount <- nCount + 1
          }
        } else {
          
          #--- if mindist is not defined ---#
          
          if (add_temp$cell %in% add_strata$cell) next
          
          add_strata <- rbind(add_strata, add_temp[, c("cell", "X", "Y", "strata", "type", "rule", extraCols)])
          
          nCount <- nCount + 1
        }
      }
    }
  }
  
  if (nCount < n) {
    message(sprintf("Strata %s: couldn't select required number of samples: %i instead of %i \n", s, nCount, n))
  }
  
  return(add_strata)
  
}

#' AHELS nSamp
#' @family rules
#' @rdname rules
#' @keywords internal
#' @export

ahels_nSamp <- function(nSamp,
                        nQuant,
                        tolerance,
                        nb,
                        underRep,
                        ratio,
                        ratOrderUnder,
                        matCovDens,
                        matCovSampDens,
                        samples,
                        mats,
                        vals){
  
  if (!is.numeric(nSamp)) {
    stop("'nSamp' must be type numeric", call. = FALSE)
  }
  
  if(tolerance > 0){
    
    message(paste0("A tolerance of ",tolerance," has been provided. Samples will be added until ",nSamp," is reached or until sampling ratios are all >= ", 1 - tolerance ,"."))
    
    
  } else {
    
    message(paste0("Samples will be added until ",nSamp," is reached or until sampling ratios are all >= 1."))
    
  }
  
  #--- begin sampling from highest discrepancy to lowest ---#
  
  newSamp <- nSamp
  sTot <- 0
  
  #--- sampling ---#
  
  while (newSamp != 0) {
    
    #--- determine the greatest discrepancy between sample and covariate data ---#
    
    repRankUnder <- which(underRep[, 3] == ratOrderUnder[1])
    
    #--- determine row and column of most under represented quantile ---#
    
    repRow <- underRep[repRankUnder, 1]
    repCol <- underRep[repRankUnder, 2]
    
    #--- if all sampling ratios in matCovSampDens are >= 1 stop adding samples ---#
    
    if(length(repRankUnder) == 0){
      
      message(paste0("Sampling ratios are all >=", 1 - tolerance,". A total of ", sTot, " samples were added."))
      
      break
    }
    
    #--- determine number of existing samples in selected quantile ---#
    sampExist <- floor(nrow(samples) * matCovSampDens[repRow, repCol])
    
    #--- determine max number of samples based on covariate density ---#
    
    sampOptim <- ceiling(nrow(samples) * matCovDens[repRow, repCol])
    
    #--- number of samples needed ---#
    
    sampNeed <- sampOptim - sampExist
    
    #--- we have a limited number of samples so we need to be sure not to over allocate ---#
    if (newSamp <= sampNeed) sampNeed <- newSamp
    
    #--- selecting covariates based on quantile chosen ---#
    
    covLower <- mats$matQ[repRow, repCol]
    
    covUpper <- mats$matQ[repRow + 1, repCol]
    
    #--- subset covariate dataset for potential new samples ---#
    
    valsSub <- vals[vals[, (2 + repCol)] >= covLower & vals[, (2 + repCol)] <= covUpper, ]
    
    #--- randomly sample within valsSub and extract randomly sampled cells ---#
    
    addSamp <- sample(nrow(valsSub), sampNeed)
    
    valsSubSamp <- valsSub[addSamp, ]
    
    valsSubSamp$type <- "new"
    
    #--- remove samples from pool to ensure same cells are not sampled again ---#
    
    vals <- vals[-addSamp, ]
    
    #--- add new samples to existing sample dataframe ---#
    
    samples <- rbind(samples, valsSubSamp)
    
    #--- update loop parameters ---#
    
    message("Quantile ", paste0("[",repRow, ",",repCol, "]"), " - A total of ", sampNeed, " samples have been allocated.")
    
    #--- update total allocated samples ---#
    
    sTot <- sTot + sampNeed
    
    #--- update available sample number ---#
    
    newSamp <- newSamp - sampNeed
    
    #--- recompute ratio's in the presence of newly added samples ---#
    
    matCovSamp <- mat_covNB(vals = samples[4:ncol(samples)], nQuant = nQuant, nb = nb, matQ = mats$matQ)
    
    #--- Change 0's to very small number to avoid division issues ---#
    
    matCovSamp[which(matCovSamp == 0)] <- 0.0000001
    
    #--- Create density matrix from covariates and length of mraster ---#
    
    matCovSampDens <- matCovSamp / nrow(samples)
    
    ### --- Selection of new samples based on density ---###
    
    #--- Ratio and ordering of data density and covariate density ---#
    
    ratio <- matCovSampDens / matCovDens
    
    #--- order the densities based on representation ---#
    
    #--- low to high ---#
    
    ratOrderUnder <- order(ratio)
    
    #--- Outline quantiles that are underrepresented (< 1) in the sample ---#
    
    underRep <- which(ratio < (1 - tolerance), arr.ind = TRUE)
    underRep <- cbind(underRep, which(ratio < (1 - tolerance)))
  }
  
  return(list(samples = samples, ratio = ratio, sTot = sTot))
  
}


#' AHELS threshold
#' @family rules
#' @rdname rules
#' @keywords internal
#' @export

ahels_threshold <- function(threshold,
                            tolerance,
                            nQuant,
                            nb,
                            ratio,
                            underRep,
                            ratOrderUnder,
                            matCovDens,
                            matCovSampDens,
                            samples,
                            mats,
                            vals){
  
  #--- Total samples added ---#
  sTot <- 0

  if(tolerance > 0){
    
    message(paste0("Threshold of ", threshold, " with a tolerance of ", tolerance ," provided. Samples will be added until sampling ratios are >= ", threshold - tolerance, "."))
    
    
  } else {
    
    message(paste0("Threshold of ", threshold, " provided. Samples will be added until sampling ratios are >= ", threshold, "."))
    
  }
  
  threshold <- threshold - tolerance

  ### --- If 'nSamp' is not provided a threshold is used ---###
  
  while (isTRUE(any(ratio < threshold))) {
    
    #--- determine the greatest discrepancy between sample and covariate data ---#
    
    repRankUnder <- which(underRep[, 3] == ratOrderUnder[1])
    
    #--- determine row and column of most under represented quantile ---#
    
    repRow <- underRep[repRankUnder, 1]
    repCol <- underRep[repRankUnder, 2]
    
    #--- determine number of existing samples in selected quantile ---#
    
    sampExist <- floor(nrow(samples) * matCovSampDens[repRow, repCol])
    
    #--- determine max number of samples based on covariate density ---#
    
    sampOptim <- ceiling(nrow(samples) * matCovDens[repRow, repCol])
    
    #--- number of samples needed ---#
    
    sampNeed <- sampOptim - sampExist
    
    #--- selecting covariates based on quantile chosen ---#
    
    covLower <- mats$matQ[repRow, repCol]
    
    covUpper <- mats$matQ[repRow + 1, repCol]
    
    #--- subset covariate dataset for potential new samples ---#
    
    valsSub <- vals[vals[, (2 + repCol)] >= covLower & vals[, (2 + repCol)] <= covUpper, ]
    
    #--- randomly sample within valsSub and extract randomly sampled cells ---#
    
    addSamp <- sample(nrow(valsSub), sampNeed)
    
    valsSubSamp <- valsSub[addSamp, ]
    
    valsSubSamp$type <- "new"
    
    #--- remove samples from pool to ensure same cells are not sampled again ---#
    
    vals <- vals[-addSamp, ]
    
    #--- add new samples to existing sample dataframe ---#
    
    samples <- rbind(samples, valsSubSamp)
    
    #--- update loop parameters ---#
    
    message("Quantile ", paste0("[",repRow, ",",repCol, "]"), " - A total of ", sampNeed, " samples have been allocated.")
    
    #--- update total allocated samples ---#
    
    sTot <- sTot + sampNeed
    
    #--- recompute ratio's in the presence of newly added samples ---#
    
    matCovSamp <- mat_covNB(vals = samples[4:ncol(samples)], nQuant = nQuant, nb = nb, matQ = mats$matQ)
    
    #--- Change 0's to very small number to avoid division issues ---#
    
    matCovSamp[which(matCovSamp == 0)] <- 0.0000001
    
    #--- Create density matrix from covariates and length of mraster ---#
    
    matCovSampDens <- matCovSamp / nrow(samples)
    
    ### --- Selection of new samples based on density ---###
    
    #--- Ratio and ordering of data density and covariate density ---#
    
    ratio <- matCovSampDens / matCovDens
    
    #--- order the densities based on representation ---#
    
    #--- low to high ---#
    
    ratOrderUnder <- order(ratio)
    
    #--- high to low ---#
    
    ratOrderOver <- rev(ratOrderUnder)
    
    #--- Outline quantiles that are underrepresented (< 1) in the sample ---#
    
    underRep <- which(ratio < threshold, arr.ind = TRUE)
    underRep <- cbind(underRep, which(ratio < threshold))
  }
  
  return(list(samples = samples, ratio = ratio, sTot = sTot))
  
}
