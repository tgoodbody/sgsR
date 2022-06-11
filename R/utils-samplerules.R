#' Stratified sampling rules
#'
#' @inheritParams sample_strat
#' @family rules
#' @name rules
NULL

#' Stratified rule 1
#' @family rules
#' @rdname rules
#' @keywords internal
#' @export

rule1 <- function(n, #number of samples
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

rule2 <- function(n, #number of samples
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

