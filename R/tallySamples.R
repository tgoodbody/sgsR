tallySamples <- function(raster,ns){
  
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
           total = as.integer(freq*ns)) %>%
    
    #--- if a value equates to <1 it will have 0 samples --- change 0 to 1 ---#
    
    mutate(total = replace(total, total == 0, 1)) %>%
    dplyr::select(strata, total) %>%
    as.data.frame()
  
  toSample

  
}