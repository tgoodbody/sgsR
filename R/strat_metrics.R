#' Stratify raster using metric quantiles
#' @family stratify functions
#'
#' @inheritParams strat_kmeans
#' @param metric Character. Name of primary metric to stratify
#' @param metric2 Character. Name of secondary metric to stratify
#' @param b Numeric. Number of desired strata for metric
#' @param b2 Numeric. Number of desired strata for metric2
#' @param samp Numeric. Determines proportion of cells to plot for strata visualization. Lower values reduce processing time.
#'
#' @return list where \code{kmeans} is all principal component analysis data and \code{raster} is the output stratification spatRaster
#' 
#' @export

strat_metrics <- function(raster,
                       metric,
                       metric2 = NULL,
                       b,
                       b2 = NULL,
                       plot = FALSE,
                       samp = 1){

  #--- error handling ---#
  if (!inherits(raster,"SpatRaster"))
    stop("all specified bands must be type SpatRaster", call. = FALSE)

  if (!is.character(metric))
    stop("'metric' must be type character")

  if (!is.numeric(b))
    stop("'b' must be type numeric")

  if (!is.logical(plot))
    stop("'plot' must be type logical")
  
  if (!is.numeric(samp))
    stop("'samp' must be type numeric")

  if (is.null(metric2)){
    if (!is.null(b2))
      message("You are stratifying with only 1 metric but specified 'b2' - ignoring.")

    #--- Extract values from raster ---#
    vals <- terra::subset(raster,metric) %>%
      terra::values()

    vals[!is.finite(vals)] <- NA

    #--- Determine index of each cell so to map values correctly without NA ---#
    idx <- !is.na(vals)

    #--- Remove NA / NaN / Inf values ---#
    df <- vals %>%
      as.data.frame() %>%
      filter(!is.na(.))

    metric <- ensym(metric)

    #--- Split metric distribution in to number specified by 'breaks' ---#
    dfc <- df %>%
      mutate(class = ntile(!!metric,b))

    #--- convert back to original raster extent ---#
    vals[idx] <- dfc$class

    #--- set newly stratified values ---#
    rout <- terra::setValues(raster[[1]],vals)
    names(rout) <- "strata"

  }

  if (!is.null(metric2)){

    if (!is.character(metric2))
      stop("'metric2' must be type character")

    if (is.null(b2))
      stop("If using 2 metrics to stratify, 'b2' must be defined")

    #--- Extract values from raster ---#
    vals <- terra::subset(raster,c(metric,metric2)) %>%
      terra::values()
    
    vals[!is.finite(vals)] <- NA

    #--- Determine index of each cell so to map values correctly without NA ---#
    idx <- is.finite(vals[,1]) & is.finite(vals[,2])

    #--- Remove NA / NaN / Inf values ---#
    df <- vals %>%
      as.data.frame() %>%
      filter(!is.na(.))

    metric <- ensym(metric)
    metric2 <- ensym(metric2)

    #--- Split metric distribution in to number specified by 'breaks' ---#
    dfc <- df %>%
      #--- define b classes ---#
      mutate(class1 = ntile(!!metric,b)) %>%
      #--- group by class to sub stratify ---#
      group_by(class1) %>%
      #--- define b2 classes ---#
      mutate(class2 = ntile(!!metric2,b2)) %>%
      #--- combine classes ---#
      group_by(class1,class2) %>%
      #--- establish newly formed unique class ---#
      mutate(class = cur_group_id())

    #--- convert back to original raster extent ---#
    vals[,1][idx] <- dfc$class

    #--- set newly stratified values ---#
    rout <- terra::setValues(raster[[1]],vals[,1])
    names(rout) <- "strata"

    if (plot == TRUE){
      if (samp > 1 | samp < 0)
        stop("'samp' must be > 0 & <= 1")

      #--- set up colour palette ---#
      ncol <- b * b2
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
                     metric = metric,
                     metric2 = metric2,
                     samp = samp)

      print(q)

    }

    return(rout)

  }

  if (plot == TRUE){

    #--- set up colour palette ---#
    ncol <- b
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

    terra::plot(rout, main = 'Classes', col=sample(col_vector, ncol))

    return(rout)

  }

}


