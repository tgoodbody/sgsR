
strat_vars <- function(raster,
                       var1,
                       var2 = NULL,
                       b1,
                       b2 = NULL,
                       plot = TRUE){

  if (!inherits(raster,"SpatRaster"))
    stop("all specified bands must be type SpatRaster", call. = FALSE)

  if (!is.character(var1))
    stop("'var1' must be type character")

  if (!is.numeric(b1))
    stop("'b1' must be type numeric")

  if (!is.logical(plot))
    stop("'plot' must be type logical")

  if (is.null(var2)){
    if (!is.null(b2))
      message("You are stratifying with only 1 variable but specified 'b2' - ignoring.")

    #--- Extract values from raster ---#
    vals <- terra::subset(raster,var1) %>%
      terra::values()

    vals[is.nan(vals)] <- NA
    vals[is.infinite(vals)] <- NA

    #--- Determine index of each cell so to map values correctly without NA ---#
    idx <- !is.na(vals)

    #--- Remove NA / NaN / Inf values ---#
    df <- vals %>%
      as.data.frame() %>%
      filter(complete.cases(.))

    var1 <- ensym(var1)

    #--- Split var1 distribution in to number specified by 'breaks' ---#
    dfc <- df %>%
      mutate(class = ntile(!!var1,b1))

    #--- convert back to original raster extent ---#
    vals[idx] <- dfc$class

    #--- set newly stratified values ---#
    rout <- terra::setValues(raster[[1]],vals)
    names(rout) <- "class"

  }

  if (!is.null(var2)){

    if (!is.character(var2))
      stop("'var2' must be type character")

    if (is.null(b2))
      stop("If using 2 variables to stratify, 'b2' must be defined")

    #--- Extract values from raster ---#
    vals <- terra::subset(raster,c(var1,var2)) %>%
      terra::values()

    #--- Determine index of each cell so to map values correctly without NA ---#
    idx <- is.finite(vals[,1]) & is.finite(vals[,2])

    #--- Remove NA / NaN / Inf values ---#
    df <- vals %>%
      as.data.frame() %>%
      filter(complete.cases(.))

    var1 <- ensym(var1)
    var2 <- ensym(var2)

    #--- Split var1 distribution in to number specified by 'breaks' ---#
    dfc <- df %>%
      #--- define b1 classes ---#
      mutate(class1 = ntile(!!var1,b1)) %>%
      #--- group by class to sub stratify ---#
      group_by(class1) %>%
      #--- define b2 classes ---#
      mutate(class2 = ntile(!!var2,b2)) %>%
      #--- combine classes ---#
      group_by(class1,class2) %>%
      #--- establish newly formed unique class ---#
      mutate(class = cur_group_id())

    #--- convert back to original raster extent ---#
    vals[,1][idx] <- dfc$class

    #--- set newly stratified values ---#
    rout <- terra::setValues(raster[[1]],vals[,1])
    names(rout) <- "class"

    if (plot == TRUE){

      #--- set up colour palette ---#
      ncol <- b1 * b2
      qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'seq',]
      col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

      terra::plot(rout, main = 'Classes', col=sample(col_vector, ncol))

      coordsgrps <- dfc %>%
        group_by(class) %>%
        arrange(class) %>%
        nest() %>%
        ungroup()

      q <- classPlot(dfc = dfc,
                     coordsgrps = coordsgrps,
                     var1 = var1,
                     var2 = var2)

      print(q)

    }

    return(rout)

  }

  if (plot == TRUE){

    #--- set up colour palette ---#
    ncol <- b1
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

    terra::plot(rout, main = 'Classes', col=sample(col_vector, ncol))

    return(rout)

  }

}


