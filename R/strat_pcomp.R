strat_pcomp <- function(raster,
                        ncp,
                        b1,
                        b2 = NULL,
                        scale = TRUE,
                        plot = TRUE,
                        samp = 1)
  
{
  #--- Error management ---#
  if (!inherits(raster,"SpatRaster"))
    stop("raster must be type SpatRaster", call. = FALSE)
  
  if (!is.numeric(ncp))
    stop("'ncp' must be type numeric")
  
  if (!is.numeric(b1))
    stop("'b1' must be type numeric")
  
  if (!is.logical(scale))
    stop("'scale' must be type logical")
  
  if (!is.logical(plot))
    stop("'plot' must be type logical")
  
  if (!is.numeric(samp))
    stop("'samp' must be type numeric")
  
  if (samp > 1.0 | samp < 0)
    stop("'samp' must be > 0 & <= 1")
  
  #--- Extract values from raster ---#
  vals <- terra::values(raster)
  
  vals[is.nan(vals)] <- NA
  vals[is.infinite(vals)] <- NA
  
  #--- Determine index of each cell so to map values correctly without NA ---#
  idx <- !is.na(vals)
  
  if ( scale == TRUE ){
    
    if ( is.null(b2) ){
      
      #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#
      pca <- FactoMineR::PCA(vals, ncp = ncp, scale.unit = TRUE,graph = FALSE)
      
      ########################################################
      ### WHICH PRINCIPAL COMPONENT METHOD SHOULD BE USED? ###
      ########################################################
      
      #--- extract PCA values ---#
      pcavals <- as.data.frame(pca$ind$coord)
      
      #--- Split PCA distribution in to number specified by 'b1' ---#
      pcagrps <- pcavals[idx,] %>%
        #--- define b1 classes ---#
        mutate(class = ntile(Dim.1,b1))
      
      #--- convert back to original raster extent ---#
      vals[idx] <- pcagrps$class
      
      #--- set newly stratified values ---#
      rout <- terra::setValues(raster[[1]],vals)
      names(rout) <- "class"
      
    } else {
      
      if (!is.numeric(b2))
        stop("'b2' must be type numeric")
      
      #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#
      pca <- FactoMineR::PCA(vals, ncp = ncp, scale.unit = TRUE,graph = FALSE)
      
      ########################################################
      ### WHICH PRINCIPAL COMPONENT METHOD SHOULD BE USED? ###
      ########################################################
      
      #--- extract PCA values ---#
      pcavals <- as.data.frame(pca$ind$coord)
      
      #--- Split PCA distribution in to number specified by 'b1' ---#
      pcagrps <- pcavals[idx,] %>%
        #--- define b1 classes ---#
        mutate(class1 = ntile(Dim.1,b1)) %>%
        #--- group by class to sub stratify ---#
        group_by(class1) %>%
        #--- define b2 classes ---#
        mutate(class2 = ntile(Dim.2,b2)) %>%
        #--- combine classes ---#
        group_by(class1,class2) %>%
        #--- establish newly formed unique class ---#
        mutate(class = cur_group_id())
      
      #--- convert back to original raster extent ---#
      vals[idx] <- pcagrps$class
      
      #--- set newly stratified values ---#
      rout <- terra::setValues(raster[[1]],vals)
      names(rout) <- "class"
      
    }
    
  } else {
    
    if ( is.null(b2) ){
      
      #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#
      pca <- FactoMineR::PCA(vals, ncp = ncp, scale.unit = FALSE,graph = FALSE)
      
      ########################################################
      ### WHICH PRINCIPAL COMPONENT METHOD SHOULD BE USED? ###
      ########################################################
      
      #--- extract PCA values ---#
      pcavals <- as.data.frame(pca$ind$coord)
      
      #--- Split PCA distribution in to number specified by 'b1' ---#
      pcagrps <- pcavals[idx,] %>%
        #--- define b1 classes ---#
        mutate(class = ntile(Dim.1,b1))
      
      #--- convert back to original raster extent ---#
      vals[idx] <- pcagrps$class
      
      #--- set newly stratified values ---#
      rout <- terra::setValues(raster[[1]],vals)
      names(rout) <- "class"
      
    } else {
      
      if (!is.numeric(b2))
        stop("'b2' must be type numeric")
      
      #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#
      pca <- FactoMineR::PCA(vals, ncp = ncp, scale.unit = FALSE,graph = FALSE)
      
      ########################################################
      ### WHICH PRINCIPAL COMPONENT METHOD SHOULD BE USED? ###
      ########################################################
      
      #--- extract PCA values ---#
      pcavals <- as.data.frame(pca$ind$coord)
      
      #--- Split PCA distribution in to number specified by 'b1' ---#
      pcagrps <- pcavals[idx,] %>%
        #--- define b1 classes ---#
        mutate(class1 = ntile(Dim.1,b1)) %>%
        #--- group by class to sub stratify ---#
        group_by(class1) %>%
        #--- define b2 classes ---#
        mutate(class2 = ntile(Dim.2,b2)) %>%
        #--- combine classes ---#
        group_by(class1,class2) %>%
        #--- establish newly formed unique class ---#
        mutate(class = cur_group_id())
      
      #--- convert back to original raster extent ---#
      vals[idx] <- pcagrps$class
      
      #--- set newly stratified values ---#
      rout <- terra::setValues(raster[[1]],vals)
      names(rout) <- "class"
      
    }
    
  }
  
  if (plot == TRUE){
    
    terra::plot(rout)
    
    coordsgrps <- pcagrps %>%
      group_by(class) %>%
      arrange(class) %>%
      na.omit() %>%
      nest() %>%
      ungroup()
    
    p <- classPlot(pcagrps,
                   coordsgrps,
                   var1 = "Dim.1",
                   var2 = "Dim.2",
                   samp)
    
    print(p)
    
  }
  
  return(rout)
  
}
