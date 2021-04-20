# raster = spatRaster. The input raster used to determine the proportional number of samples to collect based on number of strata pixels
# ns = numeric. The number of samples the user specifies to take

tallySamples <- function(raster,
                         ns,
                         existing = NULL){
  
  #--- determine crs of input raster ---#
  crs <- crs(raster)
  
  vals <- values(raster) %>% 
    as.data.frame()
  
  names(vals) <- "strata"
  
  #--- determine number of samples within each strata ---#
  toSample <- vals %>% 
    na.omit() %>%
    group_by(strata) %>% 
    summarize(n= n()) %>% 
    mutate(freq = n / sum(n),
           total = freq*ns) %>%
    
    #--- if a value equates to <1 it will have 0 samples --- change 0 to 1 ---#
    
    #########################################
    #### What other method could be used ####
    #########################################
  
    mutate(total = replace(total, total < 1, 1)) %>%
    mutate(total = round(total)) %>%
    dplyr::select(strata, total) %>%
    as.data.frame()
  
  #--- determine whether there is a difference between 'ns' and the number of samples provided based on count ---#
  
  diff <- sum(toSample$total) - ns
  
  #--- if there is a difference remove samples from classes with the most samples ---#
  
  if ( diff != 0 ) {
    
    if ( diff > 0 ) {
      
      #--- determine the largest total samples size among strata ---#
      
      maxTotal <- max(toSample$total)
      
      #--- subtract 'diff' from largest sample size ---#
      
      toSample <- toSample %>% 
        mutate(total = replace(total,
                               total == maxTotal,
                               maxTotal - abs(diff))
        )
      
    } else {
      
      #--- determine the largest total samples size among strata ---#
      
      minTotal <- min(toSample$total)
      
      #--- add 'diff' to smallest sample size ---#
      
      toSample <- toSample %>% 
        mutate(total = replace(total,
                               total == minTotal,
                               minTotal + abs(diff))
        )
      
    }
    
  }
  
  #--- if existing is provided include already sampled plots to achieve the total number ---#
  
  if (!missing(existing)) {
    
    #--- convert existing to data frame of strata values ---#
    
    existing <- data.frame(strata = existing$strata) 
    
    #--- determine number of samples for each strata ---#
    
    existing <- existing %>%
      group_by(strata) %>% 
      arrange() %>%
      summarize(eTotal= n())
    
    #--- if the strata for toSample and existing are not identical throw an error ---#
    if (!identical(unique(existing$strata),unique(toSample$strata)))
      stop("Strata for 'raster' and 'existing' are not identical")
    
    #--- join the 2 df together and subtrace the number of existing plots by strata from toSample ---#
    toSample <- toSample %>%
      left_join(existing, by = "strata") %>%
      mutate(total = total - eTotal) %>%
      dplyr::select(-eTotal) %>%
      as.data.frame()

  }
  
  toSample

}
