# raster = spatRaster. Multiband ALS metrics raster.
# ncp = Character. Number of components to create.
# b1 = Numeric. Number of desired strata for first principal component.
# b2 = Numeric. Number of desired strata for second principal component.
# scale - Logical. Determines whether centering and scaling of data should be conducted prior to principal component analysis.
# plot = Logical. Plots output strata raster and visualized strata with boundary dividers.
# samp = Numeric. Determines proportion of cells to plot for strata visualization. Lower values reduce processing time.

## output is a list where:
#'$pca' is all principal component analysis data
#'$raster' is the output stratification spatRaster

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
  
  #--- Extract values from raster ---#
  vals <- terra::values(raster)
  
  vals[!is.finite(vals)] <- NA
  
  #--- Determine index of each cell so to map values correctly without NA ---#
  idx <- !is.na(vals)
  
  if ( scale == TRUE ){
    
    if ( is.null(b2) ){
      
      #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#
      pca <- suppressWarnings(FactoMineR::PCA(vals, ncp = ncp, scale.unit = TRUE,graph = FALSE))
      
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
      names(rout) <- "strata"
      
    } else {
      
      if (!is.numeric(b2))
        stop("'b2' must be type numeric")
      
      #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#
      pca <- suppressWarnings(FactoMineR::PCA(vals, ncp = ncp, scale.unit = TRUE,graph = FALSE))
      
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
      names(rout) <- "strata"
      
    }
    
  } else {
    
    if ( is.null(b2) ){
      
      #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#
      pca <- suppressWarnings(FactoMineR::PCA(vals, ncp = ncp, scale.unit = FALSE,graph = FALSE))
      
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
      names(rout) <- "strata"
      
    } else {
      
      if (!is.numeric(b2))
        stop("'b2' must be type numeric")
      
      #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#
      pca <- suppressWarnings(FactoMineR::PCA(vals, ncp = ncp, scale.unit = FALSE,graph = FALSE))
      
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
      names(rout) <- "strata"
      
    }
    
  }
  
  if (plot == TRUE){
    if (samp > 1 | samp < 0)
      stop("'samp' must be > 0 & <= 1")
    
    terra::plot(rout)
    
    coordsgrps <- pcagrps %>%
      group_by(class) %>%
      arrange(class) %>%
      na.omit() %>%
      nest() %>%
      ungroup()
    
    p <- classPlot(pcagrps,
                   coordsgrps,
                   metric = "Dim.1",
                   metric2 = "Dim.2",
                   samp)
    
    print(p)
    
  }
  
 out <- list(pca = pca, raster = rout)
 
 out
  
}
