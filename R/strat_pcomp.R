#' Stratify metrics raster using principal components and quantile breaks
#' @family stratify functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams strat_quantiles
#' 
#' @param ...
#' 
#' @importFrom magrittr %>%
#' @importFrom methods is
#' 
#' @return output stratification \code{spatRaster}, or a list when \code{details = TRUE}.
#' 
#' @export

strat_pcomp <- function(mraster,
                        nstrata,
                        nstrata2 = NULL,
                        scale = TRUE,
                        plot = FALSE,
                        samp = 1,
                        details = FALSE,
                        ...)
  
{
  
  #--- set global vars ---#
  
  raster <- class1 <- class2 <- Dim.1 <- Dim.2 <- NULL
  
  #--- Error management ---#
  
  if (!inherits(mraster,"SpatRaster"))
    stop("mraster must be type SpatRaster", call. = FALSE)
  
  if (!is.numeric(nstrata))
    stop("'nstrata' must be type numeric")
  
  if (!is.logical(scale))
    stop("'scale' must be type logical")
  
  if (!is.logical(plot))
    stop("'plot' must be type logical")
  
  if (!is.numeric(samp))
    stop("'samp' must be type numeric")
  
  if (!is.logical(details))
    stop("'details' must be type logical")
  
  #--- Extract values from mraster ---#
  
  vals <- terra::values(mraster)
  
  vals[!is.finite(vals)] <- NA
  
  #--- Determine index of each cell so to map values correctly without NA ---#
  
  idx <- !is.na(vals)
  
  if ( scale == TRUE ){
    
    if ( is.null(nstrata2) ){
      
      #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#
      
      pca <- suppressWarnings(FactoMineR::PCA(vals, ncp = 2, scale.unit = TRUE, graph = FALSE, ...))
      
      ########################################################
      ### WHICH PRINCIPAL COMPONENT METHOD SHOULD BE USED? ###
      ########################################################
      
      #--- extract PCA values ---#
      
      pcavals <- as.data.frame(pca$ind$coord)
      
      #--- Split PCA distribution into number specified by 'nstrata' ---#
      
      pcagrps <- pcavals[idx,] %>%
        #--- define nstrata classes ---#
        dplyr::mutate(class = dplyr::ntile(Dim.1,nstrata))
      
      #--- convert back to original mraster extent ---#
      
      vals[idx] <- pcagrps$class
      
      #--- set newly stratified values ---#
      
      rout <- terra::setValues(mraster[[1]],vals)
      names(rout) <- "strata"
      
    } else {
      
      if (!is.numeric(nstrata2))
        stop("'nstrata2' must be type numeric")
      
      #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#
      
      pca <- suppressWarnings(FactoMineR::PCA(vals, ncp = 2, scale.unit = TRUE, graph = FALSE, ...))
      
      ########################################################
      ### WHICH PRINCIPAL COMPONENT METHOD SHOULD BE USED? ###
      ########################################################
      
      #--- extract PCA values ---#
      
      pcavals <- as.data.frame(pca$ind$coord)
      
      #--- Split PCA distribution in to number specified by 'nstrata' ---#
      
      pcagrps <- pcavals[idx,] %>%
        #--- define nstrata classes ---#
        dplyr::mutate(class1 = dplyr::ntile(Dim.1,nstrata)) %>%
        #--- group by class to sub stratify ---#
        dplyr::group_by(class1) %>%
        #--- define nstrata2 classes ---#
        dplyr::mutate(class2 = dplyr::ntile(Dim.2,nstrata2)) %>%
        #--- combine classes ---#
        dplyr::group_by(class1,class2) %>%
        #--- establish newly formed unique class ---#
        dplyr::mutate(class = dplyr::cur_group_id())
      
      #--- convert back to original raster extent ---#
      
      vals[idx] <- pcagrps$class
      
      #--- set newly stratified values ---#
      
      rout <- terra::setValues(mraster[[1]],vals)
      names(rout) <- "strata"
      
    }
    
  } else {
    
    if ( is.null(nstrata2) ){
      
      #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#
      
      pca <- suppressWarnings(FactoMineR::PCA(vals, ncp = 2, scale.unit = FALSE, graph = FALSE, ...))
      
      ########################################################
      ### WHICH PRINCIPAL COMPONENT METHOD SHOULD BE USED? ###
      ########################################################
      
      #--- extract PCA values ---#
      
      pcavals <- as.data.frame(pca$ind$coord)
      
      #--- Split PCA distribution in to number specified by 'nstrata' ---#
      
      pcagrps <- pcavals[idx,] %>%
        #--- define nstrata classes ---#
        dplyr::mutate(class = dplyr::ntile(Dim.1,nstrata))
      
      #--- convert back to original raster extent ---#
      
      vals[idx] <- pcagrps$class
      
      #--- set newly stratified values ---#
      
      rout <- terra::setValues(raster[[1]],vals)
      names(rout) <- "strata"
      
    } else {
      
      if (!is.numeric(nstrata2))
        stop("'nstrata2' must be type numeric")
      
      #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#
      
      pca <- suppressWarnings(FactoMineR::PCA(vals, ncp = 2, scale.unit = FALSE, graph = FALSE, ...))
      
      ########################################################
      ### WHICH PRINCIPAL COMPONENT METHOD SHOULD BE USED? ###
      ########################################################
      
      #--- extract PCA values ---#
      
      pcavals <- as.data.frame(pca$ind$coord)
      
      #--- Split PCA distribution in to number specified by 'nstrata' ---#
      
      pcagrps <- pcavals[idx,] %>%
        #--- define nstrata classes ---#
        dplyr::mutate(class1 = dplyr::ntile(Dim.1,nstrata)) %>%
        #--- group by class to sub stratify ---#
        dplyr::group_by(class1) %>%
        #--- define b2 classes ---#
        dplyr::mutate(class2 = dplyr::ntile(Dim.2,nstrata2)) %>%
        #--- combine classes ---#
        dplyr::group_by(class1,class2) %>%
        #--- establish newly formed unique class ---#
        dplyr::mutate(class = dplyr::cur_group_id())
      
      #--- convert back to original raster extent ---#
      
      vals[idx] <- pcagrps$class
      
      #--- set newly stratified values ---#
      
      rout <- terra::setValues(mraster[[1]],vals)
      names(rout) <- "strata"
      
    }
    
  }
  
  if (isTRUE(plot)){
    if (samp > 1 | samp < 0)
      stop("'samp' must be > 0 & <= 1")
    
    terra::plot(rout)
    
    coordsgrps <- pcagrps %>%
      dplyr::group_by(class) %>%
      dplyr::arrange(class) %>%
      stats::na.omit() %>%
      tidyr::nest() %>%
      dplyr::ungroup()
    
    p <- classPlot(pcagrps,
                   coordsgrps,
                   metric = "Dim.1",
                   metric2 = "Dim.2",
                   samp)
    
    print(p)
    
  }
  
  #--- Output based on 'details' to return raster alone or list with details ---#
  
  if ( isTRUE(details) ){
    
    #--- create list to assign pca info and output raster ---#
    
    out <- list(details = pca, raster = rout)
    
    return(out)
    
    
  } else {
    
    #--- just output raster ---#
    
    return(rout)
    
  }

  
}
