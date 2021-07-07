#' Calculate and rasterize principal components from a metric raster
#' 
#' @family calculate functions
#' 
#' @inheritParams sample_srs
#' @inheritParams strat_kmeans
#' 
#' @param nComp Numeric. Value indicating number of principal components to be rasterized.
#' 
#'
#' @importFrom magrittr %>%
#' @importFrom methods is
#' 
#' @return Output raster with specified number of principal components as layers.
#' 
#' @export

calculate_pcomp <- function(mraster = NULL,
                            nComp = NULL,
                            scale = TRUE,
                            plot = FALSE,
                            details = FALSE,
                            ...)
{
  
  #--- error handling ---#
  
  if (!inherits(mraster, "SpatRaster"))
    stop("'raster' must be type SpatRaster", call. = FALSE)
  
  if (!is.numeric(nComp))
    stop("'ncomp' must be type numeric")
  
  if (!is.logical(scale))
    stop("'scale' must be type logical")
  
  if (!is.logical(plot))
    stop("'plot' must be type logical")
  
  if (nComp > terra::nlyr(mraster))
    stop("nComp must be <= the number of layers in 'mraster'")
  
  #--- Extract values from mraster ---#
  
  vals <- terra::values(mraster)
  
  vals[!is.finite(vals)] <- NA
  
  #--- Determine index of each cell so to map values correctly without NA ---#
  
  idx <- !is.na(vals)
  
  #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#
  
  PCA <- suppressWarnings(FactoMineR::PCA(vals, ncp = nComp, scale.unit = scale, graph = FALSE))
  
  ########################################################
  ### WHICH PRINCIPAL COMPONENT METHOD SHOULD BE USED? ###
  ########################################################

  #--- extract cell level pca values ---#
  pcavals <- as.data.frame(PCA$ind$coord)
  
  #--- create loop to allocate values to cell level taking into account potential NA using 'idx' ---#
  rs <- list()
  
  for(i in 1:nComp){
    
    vals[idx] <- pcavals[idx,i]
    
    rs[[i]] <- terra::setValues(mraster[[1]], vals)
    
  }
  
  #--- stack pca rasters ---#
  
  pcaRout <- rast(rs)
  
  #--- rename ---#
  
  names(pcaRout) <- rep(paste0("PC",seq(1,nComp,1)))
  
  #--- plot scree plot ---#
  
  if( isTRUE(plot) ){
    
    #--- visualize scree plot ---#
    print(factoextra::fviz_screeplot(PCA))
    
  }
  
  if ( isTRUE(details) ){
    
    #--- create list to assign pca info and output raster ---#
    
    out <- list(pca = PCA, raster = pcaRout)
    
    return(out)
    
    
  } else {
    
    #--- just output raster ---#
    
    return(pcaRout)
    
  }
  
}

