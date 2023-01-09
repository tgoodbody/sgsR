#' coobs algorithm sampling
#'
#' @description Perform the COunt of OBServations (coobs) algorithm using existing site data
#' and raster metrics. This algorithm aids the user in determining where additional samples
#' could be located by comparing existing samples to each pixel and associated covariates.
#' The output coobs raster could be used to constrain clhs sampling to areas that are underreprented.
#' @family calculate functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams extract_strata
#' @inheritParams strat_breaks
#'
#' @param mraster spatRaster. ALS metrics raster. Requires at least 2 layers to calculate covariance matrix.
#' @param threshold Numeric. Proxy maximum pixel quantile to avoid outliers. \code{default = 0.95}.
#' @param cores Numeric. Number of cores to use for parallel processing. \code{default = 1}.
#'
#' @references
#' Malone BP, Minansy B, Brungard C. 2019. Some methods to improve the utility of conditioned Latin hypercube sampling. PeerJ 7:e6451 DOI 10.7717/peerj.6451
#'
#' @return Output raster with coobs and classified coobs layers.
#'
#' @examples
#' \dontrun{
#' #--- Load raster and existing plots---#
#' r <- system.file("extdata", "mraster.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' e <- system.file("extdata", "existing.shp", package = "sgsR")
#' e <- sf::st_read(e)
#'
#' calculate_coobs(
#'   mraster = mr,
#'   existing = e,
#'   cores = 4
#' )
#' }
#' @note
#' Special thanks to Dr. Brendan Malone for the original implementation of this algorithm.
#'
#' @author Tristan R.H. Goodbody
#'
#' @export

calculate_coobs <- function(mraster,
                            existing,
                            cores = 1,
                            threshold = 0.95,
                            plot = FALSE,
                            filename = NULL,
                            overwrite = FALSE) {
  
  #--- check for required packages ---#
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("Package \"doParallel\" is needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  
  if (!requireNamespace("doSNOW", quietly = TRUE)) {
    stop("Packag \"doSNOW\" is needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("Package \"foreach\" is needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  
  if (!requireNamespace("snow", quietly = TRUE)) {
    stop("Package \"snow\" is needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  
  if (!requireNamespace("Rfast", quietly = TRUE)) {
    stop("Package \"Rfast\" is needed for this function to work. Please install it.",
         call. = FALSE
    )
  }
  
  #--- set global vars ---#
  
  i <- geometry <- NULL
  
  #--- Error handling ---#
  
  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster")
  }
  
  if (!inherits(existing, "data.frame") && !inherits(existing, "sf")) {
    stop("'existing' must be a data.frame or sf object")
  }
  
  if (!is.numeric(threshold)) {
    stop("'threshold' must be type numeric")
  }
  
  nb <- terra::nlyr(mraster)
  
  if (nb < 2) {
    stop("'mraster' only has 1 band. Need at least 2 bands to calculate covariance matrix")
  }
  
  #--- extract covariates data from mraster ---#
  
  vals <- terra::as.data.frame(mraster, xy = TRUE, row.names = FALSE)
  
  #--- Remove NA / NaN / Inf values ---#
  
  vals <- vals %>%
    stats::na.omit()
  
  #--- Generate covariance matrix ---#
  
  covMat <- as.matrix(stats::cov(vals[, 3:ncol(vals)]))
  
  #--- select geometry attribute ---#
  
  existing <- existing %>%
    dplyr::select(geometry)
  
  #--- extract covariates at existing sample locations ---#
  
  samples <- sgsR::extract_metrics(mraster, existing, data.frame = TRUE)
  
  #--- check if samples fall in areas where stratum values are NA ---#
  
  if(any(!complete.cases(samples))){
    
    samples_NA <- samples %>%
      dplyr::filter(!complete.cases(.)) %>%
      dplyr::mutate(type = "existing")

    samples <- samples %>%
      stats::na.omit()
    
  }
  
  #--- create parallel processing structure ---#
  
  cl <- snow::makeCluster(spec = cores, type = "SOCK")
  doSNOW::registerDoSNOW(cl)
  
  #--- create progress text bar ---#
  
  iterations <- nrow(vals)
  pb <- utils::txtProgressBar(max = iterations, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  M <- as.matrix(vals[, 3:ncol(vals)])
  sample <- nrow(M) > 10000
  if (sample) ix <- sample(1:nrow(M), 10000)
  
  `%dopar%` <- foreach::`%dopar%`
  
  #--- iterate parallel processing of mahalanobis distance ---#
  
  loop <- foreach::foreach(i = 1:iterations, .combine = "c", .options.snow = opts) %dopar% {
    cell <- vals[i, 3:ncol(vals)]
    
    #--- Determine distance for each pixel in raster ---#
    
    pixDist <- Rfast::mahala(x = M, mu = as.matrix(cell), s = covMat)
    
    
    #--- Determine min and max distance values for each ---#
    
    pixMin <- min(pixDist)
    if (sample)
      pixMax <- stats::quantile(pixDist[ix], probs = threshold)
    else
      pixMax <- stats::quantile(pixDist, probs = threshold)
    
    #--- Determine distance for each sample location ---#
    
    sampDist <- Rfast::mahala(x = as.matrix(samples[, 3:ncol(samples)]), mu = as.matrix(cell), s = covMat) # calculate distance of observations to all other pixels
    
    #--- Normalize distance between data and samples ---#
    
    sampNDist <- (sampDist - pixMin) / (pixMax - pixMin)
    
    #--- If sampDist > 1 sampDist > maxDist ---#
    
    sampNDist[sampNDist > 1] <- 1
    
    #--- larger values equate to more similarity ---#
    
    sampNDist <- 1 - sampNDist
    
    #--- establish count above threshold ---#
    
    sum(sampNDist >= threshold)
  }
  
  close(pb)
  
  #--- End parallel ---#
  snow::stopCluster(cl)
  
  #--- Coerce output from parallel to a new attribute in covariates ---#
  
  vals$nSamp <- loop
  
  #--- convert nSamp to raster ---#
  
  r <- terra::rast(as.matrix(vals[, c("x", "y", "nSamp")]), type = "xyz")
  names(r) <- "COOB"
  
  #--- classify raster into breaks on the fly ---#
  
  breaks <- unique(floor(seq(min(vals$nSamp), max(vals$nSamp), length.out = 8)))
  
  rc <- terra::classify(r, breaks, include.lowest = TRUE, right = FALSE)
  names(rc) <- "COOBclass"
  
  #--- stack 2 rasters for output ---#
  
  rout <- c(r, rc)
  
  #--- Plot output ---#
  
  if (isTRUE(plot)) {
    
    #--- plot ---#
    
    terra::plot(rc)
    terra::plot(existing, add = TRUE)
  }
  
  if (!is.null(filename)) {
    terra::writeRaster(x = rout, filename = filename, overwrite = overwrite)
    message("Output raster written to disc.")
  }
  
  return(rout)
}
