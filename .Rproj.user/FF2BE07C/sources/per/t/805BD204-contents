stat_pcomp <- function(raster,
                       nComp,
                       b1,
                       b2 = NULL,
                       spca = TRUE,
                       plot = TRUE)

  {
  #--- Error management ---#
  if (!inherits(raster,"SpatRaster"))
    stop("raster must be type SpatRaster", call. = FALSE)

  if (!is.numeric(ncomp))
    stop("'ncomp' must be type numeric")

  if (!is.numeric(b1))
    stop("'b1' must be type numeric")

  if (!is.numeric(b2))
    stop("'b2' must be type numeric")

  if (!is.logical(spca))
    stop("'spca' must be type logical")

  if (!is.logical(plot))
    stop("'plot' must be type logical")

  if ( spca == TRUE){

    if (is.null(b2)){

      #--- perform PCA using rasterPCA -- requires conversion to raster* format ---#
      pca <- RStoolbox::rasterPCA(as(raster,"Raster"), nComp = nComp, spca = TRUE)

      #--- extract PCA values ---#
      pcavals <- data.frame(values(pca$map))

      #--- Split PCA distribution in to number specified by 'b1' ---#
      pcagrps <- pcavals %>%
        #--- define b1 classes ---#
        mutate(class = ntile(PC1,b1)) %>%


      #--- convert back to original raster extent ---#
      pcavals[,1] <- pcagrps$class

      #--- set newly stratified values ---#
      rout <- terra::setValues(raster[[1]],pcavals[,1])
      names(rout) <- "class"

    } else {

      pca <- RStoolbox::rasterPCA(as(raster,"Raster"), nComp = nComp, spca = TRUE)

      pcavals <- data.frame(values(pca$map))

      #--- Split PCA1 distribution in to number specified by 'breaks' ---#
      pcagrps <- pcavals %>%
        #--- define b1 classes ---#
        mutate(class1 = ntile(PC1,b1)) %>%
        #--- group by class to sub stratify ---#
        group_by(class1) %>%
        #--- define b2 classes ---#
        mutate(class2 = ntile(PC2,b2)) %>%
        #--- combine classes ---#
        group_by(class1,class2) %>%
        #--- establish newly formed unique class ---#
        mutate(class = cur_group_id())

      #--- convert back to original raster extent ---#
      pcavals[,1] <- pcagrps$class

      #--- set newly stratified values ---#
      rout <- terra::setValues(raster[[1]],pcavals[,1])
      names(rout) <- "class"

    }

  } else {

    if (is.null(b2)){

      pca <- RStoolbox::rasterPCA(as(raster,"Raster"), nComp = nComp, spca = FALSE)

      pcavals <- data.frame(values(pca$map))

      #--- Split PCA1 distribution in to number specified by 'breaks' ---#
      pcagrps <- pcavals %>%
        #--- define b1 classes ---#
        mutate(class = ntile(PC1,b1)) %>%


        #--- convert back to original raster extent ---#
        pcavals[,1] <- pcagrps$class

      #--- set newly stratified values ---#
      rout <- terra::setValues(raster[[1]],pcavals[,1])
      names(rout) <- "class"

    } else {

      pca <- RStoolbox::rasterPCA(as(raster,"Raster"), nComp = nComp, spca = FALSE)

      pcavals <- data.frame(values(pca$map))

      #--- Split PCA1 distribution in to number specified by 'breaks' ---#
      pcagrps <- pcavals %>%
        #--- define b1 classes ---#
        mutate(class1 = ntile(PC1,b1)) %>%
        #--- group by class to sub stratify ---#
        group_by(class1) %>%
        #--- define b2 classes ---#
        mutate(class2 = ntile(PC2,b2)) %>%
        #--- combine classes ---#
        group_by(class1,class2) %>%
        #--- establish newly formed unique class ---#
        mutate(class = cur_group_id())

      #--- convert back to original raster extent ---#
      pcavals[,1] <- pcagrps$class

      #--- set newly stratified values ---#
      rout <- terra::setValues(raster[[1]],pcavals[,1])
      names(rout) <- "class"

    }

  }

  if (plot == TRUE){

    plot(pca$map)

    coordsgrps <- pcagrps %>%
      group_by(class) %>%
      arrange(class) %>%
      nest() %>%
      ungroup()

    p <- classPlot(pcagrps,
                   coordsgrps,
                   var1 = "PC1",
                   var2 = "PC2")

    print(p)

  }

  ##



}
