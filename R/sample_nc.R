#' Nearest centroid (NC) sampling
#'
#' @description Sampling using the nearest centroid (NC) approach described in Melville & Stone (2016).
#'
#' @family sample functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams sample_srs
#' 
#' @param iter Numeric. The maximum number of kmeans iterations allowed.
#' @param algorithm Character. \code{Lloyd} (default) or
#' \code{MacQueen} kmeans algorithms.
#' @param k Numeric. The number of nearest neighbours to take for each k-means center.
#' When \code{k = 1} (default), the output number of samples will match \code{nSamp}. 
#' Increases to \code{k} results in a multiplicative result total number of samples \code{nSamp * k}.
#'
#' @return An sf object with \code{nSamp} randomly sampled points.
#' 
#' @note
#' When \code{details = TRUE}, a list is returned where:
#' \enumerate{
#' \item \code{samples} output nearest centroid samples with `kcenter` attribute linking 
#' to associated kmeans centers.
#' \item \code{kmeans} is a list output of the \code{\link[stats]{kmeans}} function
#' \item \code{centers} Un-scaled kmeans center values for each layer in \code{mraster} 
#' with `kcenter` attribute to link with the same attribute in \code{samples}.
#' \item \code{kplot} is a \code{ggplot} scatterplot object visualizing the kmeans centers
#'  and associated nearest neighbor samples.
#' }
#'
#' @examples
#' #--- Load raster and access files ---#
#' r <- system.file("extdata", "mraster.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' #--- perform simple random sampling ---#
#' sample_nc(
#'   mraster = mr,
#'   nSamp = 5,
#' )
#'
#' @author Tristan R.H. Goodbody
#' 
#' @references
#' G. Melville & C. Stone (2016) Optimising nearest neighbour information—a simple,
#' efficient sampling strategy for forestry plot imputation using remotely sensed data,
#' Australian Forestry, 79:3, 217-228, DOI: 10.1080/00049158.2016.1218265
#'
#' @export

sample_nc <- function(mraster,
                      nSamp,
                      k = 1,
                      iter = 500,
                      algorithm = "Lloyd",
                      access = NULL,
                      buff_inner = NULL,
                      buff_outer = NULL,
                      plot = FALSE,
                      details = FALSE,
                      filename = NULL,
                      overwrite = FALSE){
  
  #--- set global vars ---#
  
  type <- NULL
  
  #--- check for required packages ---#
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package \"RANN\" is needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  
  #--- error handling ---#
  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster.", call. = FALSE)
  }
  
  if (!is.numeric(nSamp)) {
    stop("'nSamp' must be type numeric.", call. = FALSE)
  }
  
  if (!is.numeric(k)) {
    stop("'k' must be type numeric.", call. = FALSE)
  }
  
  if (!is.numeric(iter)) {
    stop("'iter' must be type numeric.", call. = FALSE)
  }
  
  if (!is.character(algorithm)) {
    stop("'algorithm' must be type character.", call. = FALSE)
  }
  
  if (algorithm != "Lloyd" && algorithm != "MacQueen") {
    stop("Unknown algorithm '", algorithm, "' selected. Please use 'Lloyd' (default) or 'MacQueen'.", call. = FALSE)
  }
  
  if (!is.logical(plot)) {
    stop("'plot' must be type logical.", call. = FALSE)
  }

  if (!is.logical(details)) {
    stop("'details' must be type logical.", call. = FALSE)
  }
  
  mrasterP <- mraster[[1]]
  
  #--- determine crs of input raster ---#
  crs <- terra::crs(mraster, proj = TRUE)
  
  #--- Determine number of layers in mraster ---#
  
  nc <- terra::nlyr(mraster)
  
  #--- constrain sampling by access if desired ---#
  if (!is.null(access)) {
  
    access_buff <- mask_access(raster = mraster, access = access, buff_inner = buff_inner, buff_outer = buff_outer)
    
    mraster <- access_buff$rast
  }
  
  #--- Extract values from mraster ---#
  
  vals <- terra::as.data.frame(mraster, xy = TRUE) %>%
    na.omit()

  #--- Determine index of each cell so to map values correctly without NA ---#

  d <- within(vals, valscs <- scale(vals[,3:ncol(vals)]))
  
  #--- conduct unsupervised k-means with center/scale parameters based on algorithm ---#
  
  message("K-means being performed on ", terra::nlyr(mraster), " layers with ", nSamp, " centers.")
  
  km_clust <- stats::kmeans(d[,3+nc], centers = nSamp, iter.max = iter, algorithm = algorithm)
  
  #--- perform nearest neighbour search for each kmeans center ---#
  
  idx <- RANN::nn2(data = d$valscs, query = km_clust$centers, k = k)
  
  #--- assign the 'kcenter' attribute to know which sample is the nn for which center in km_clust ---# 
  
  samples <- do.call(rbind,apply(idx$nn.idx, MARGIN = 2, FUN = function(x) as.data.frame(d[x,])[,1:(2+nc)]))
  
  samples$kcenter <- 1:nSamp
  
  #--- convert coordinates to a spatial points object ---#
  samples <- samples %>%
    sf::st_as_sf(., coords = c("x", "y"))
  
  #--- assign raster crs to spatial points object ---#
  sf::st_crs(samples) <- crs
  
  #--- rescale kmeans centers for details and plotting ---#
  
  if(nc != 1){
    
    centers <- as.data.frame(t(apply(km_clust$centers, 1, function(r)r*attr(d$valscs,'scaled:scale') + attr(d$valscs, 'scaled:center'))))
    
    
  } else {
    
    centers <- as.data.frame(km_clust$centers * attr(d$valscs,'scaled:scale') + attr(d$valscs, 'scaled:center'))
    
  }
  #--- apply names to k means clkusters centers and assign center number to match 'samples' ---#
  names(centers) <- names(mraster)
  
  centers$kcenter <- 1:nSamp
  
  #--- plotting ---#
  if(isTRUE(plot)){
    if (!is.null(access)) {
      
      #--- plot input raster and random samples ---#
      terra::plot(mrasterP[[1]])
      suppressWarnings(terra::plot(access_buff$buff, add = T, border = c("gray30"), col = "gray10", alpha = 0.1))
      suppressWarnings(terra::plot(samples, add = T, col = "black"))
    } else {
      
      #--- plot input raster and random samples ---#
      terra::plot(mrasterP[[1]])
      suppressWarnings(terra::plot(samples, add = T, col = "black"))
    }
  }
  
  #--- write samples to disk ---#
  if (!is.null(filename)) {
    
    if (!is.character(filename)) {
      stop("'filename' must be a file path character string.", call. = FALSE)
    }
    
    if (!is.logical(overwrite)) {
      stop("'overwrite' must be type logical.", call. = FALSE)
    }
    
    if (file.exists(filename) & isFALSE(overwrite)) {
      stop(paste0("'",filename, "' already exists and overwrite = FALSE."), call. = FALSE)
    }
    
    sf::st_write(samples, filename, delete_layer = overwrite)
    message("Output samples written to disc.")
  }
  
  #--- output sampling details ---#
  if (isTRUE(details)) {
    
    #--- generate basic plot showing centers and samples cells ---#
    
    #--- if raster has more than 2 layers to plot ---#
    if(nc > 1){
      
      #--- take random sample of 1000 pixels for visualization ---#
      if(nrow(vals) > 1000) vals <- vals %>% dplyr::slice_sample(n=1000)
      
      cs <- rbind(
        centers %>% 
          dplyr::mutate(type = "centers"), 
        samples %>% 
          sf::st_drop_geometry() %>%
          dplyr::mutate(type = "samples")
        )
      
      nm <- as.character(names(vals))
      
      metric <- nm[3]
      metric2 <- nm[4]
      
      metric <- ggplot2::ensym(metric)
      metric2 <- ggplot2::ensym(metric2)
      
      p <- ggplot2::ggplot(data = vals, ggplot2::aes(x = !!metric, y = !!metric2)) + 
        ggplot2::geom_point() + 
        ggplot2::geom_point(data = cs %>% dplyr::filter(type == "centers"), ggplot2::aes(x = !!metric, y = !!metric2), colour = "#F8766D", size = 2.4) +
        ggplot2::geom_point(data = cs %>% dplyr::filter(type == "samples"), ggplot2::aes(x = !!metric, y = !!metric2), colour = "#619CFF") +
        ggtitle("Nearest centroid (NC) sampling", subtitle = "Red points represent kmeans centroids, blue point represent NN samples.")

    }
    
    #--- create list to assign samples, k-means info, un-scaled cluster means, and visualization plot ---#
    
    out <- list(samples = samples, 
                kmeans = km_clust, 
                centers = centers, 
                kplot = if (exists("p")) p)
    
    return(out)
    
  } else {
    
    #--- just output raster ---#
    
    return(samples)
  }
  
}
