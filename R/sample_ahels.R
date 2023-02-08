#' Adapted Hypercube Evaluation of a Legacy Sample (ahels)
#'
#' @description Perform the adapted Hypercube Evaluation of a Legacy Sample (ahels) algorithm using
#' existing site data and raster metrics. New samples are allocated based on quantile ratios between
#' the existing sample and covariate dataset.
#'
#' @family sample functions
#'
#' @inheritParams sample_systematic
#' @inheritParams strat_kmeans
#' @inheritParams extract_strata
#'
#' @param nQuant Numeric. Number of quantiles to divide covariates and samples into. Quantiles that do not
#' cover at least 1 percent of the area of interest will be excluded and be returned as \code{NA}.
#' @param nSamp Numeric. Maximum number of new samples to allocate.
#' @param threshold Numeric. Sample quantile ratio threshold. After the threshold \code{default = 0.9} is reached,
#' no additional samples will be added. Values close to 1 can cause the algorithm to continually loop.
#' @param tolerance Numeric. Allowable tolerance (<= 0.1 (10%)) around quantile density of 1. If \code{nSamp} is used samples will be
#' added until the \code{1 - tolerance} density is reached. If \code{threshold} is used, samples will be added until the 
#' \code{threshold - tolerance} value is reached. This parameter allows the user to define a buffer around desired quantile densities
#' to permit the algorithm to not add additional samples if quantile density is very close to 1, or user-defined \code{threshold}.
#' @param matrices List. Quantile and covariance matrices generated from \code{calculate_pop(mraster = mraster, nQuant = nQuant)}.
#' Both \code{mraster} & \code{nQuant} inputs must be the same to supply the covariance matrix. Supplying the matrix allows users
#' with very large \code{mrasters} to pre-process the covariance matrix to avoid longer sampling processing times.
#' @param plot Logical. Plots samples of type \code{existing} (if provided; crosses) and \code{new} (circles) along with \code{mraster}.
#'
#' @references
#' Malone BP, Minansy B, Brungard C. 2019. Some methods to improve the utility of conditioned Latin hypercube sampling. PeerJ 7:e6451 DOI 10.7717/peerj.6451
#'
#' @return Returns sf point object with existing samples and supplemental samples added by the ahels algorithm.
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
#' sample_ahels(
#'   mraster = mr,
#'   existing = e,
#'   plot = TRUE
#' )
#'
#' #--- supply quantile and covariance matrices ---#
#' mat <- calculate_pop(mraster = mr)
#' 
#' sample_ahels(
#'   mraster = mr,
#'   existing = e,
#'   matrices = mat,
#'   nSamp = 300
#' )
#' }
#' @note
#'
#' Messages in the algorithm will state that samples have been added to under-represented quantiles. The number between
#' square brackets that follow represent the matrix row and column respectively that can be printed using \code{details = TRUE}.
#' 
#' In some cases, generally when a single metric is used as \code{mraster}, sampling ratios all be >= 1 before the 
#' \code{nSamp} number of samples are allocated. The algorithm will stop in this scenario.
#'
#' Special thanks to Dr. Brendan Malone for the original implementation of this algorithm.
#'
#' @author Tristan R.H. Goodbody
#'
#' @export

sample_ahels <- function(mraster,
                         existing,
                         nQuant = 10,
                         nSamp = NULL,
                         threshold = 0.9,
                         tolerance = 0,
                         matrices = NULL,
                         plot = FALSE,
                         details = FALSE,
                         filename = NULL,
                         overwrite = FALSE) {

  #--- Set global vars ---#

  x <- y <- X <- Y <- n <- type <- geometry <- extraCols <- NULL

  #--- Error handling ---#

  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster.", call. = FALSE)
  }

  if (!inherits(existing, "data.frame") && !inherits(existing, "sf")) {
    stop("'existing' must be a data.frame or sf object.", call. = FALSE)
  }

  if (!is.numeric(nQuant)) {
    stop("'nQuant' must be type numeric.", call. = FALSE)
  }

  if (!is.numeric(threshold)) {
    stop("'threshold' must be type numeric.", call. = FALSE)
  }
  
  if (threshold < 0 | threshold > 1) {
    stop("'threshold' must be > 0 and < 1.", call. = FALSE)
  }
  
  if (!is.logical(details)) {
    stop("'details' must be type logical.", call. = FALSE)
  }
  
  #--- tolerance ---#
  if (!is.numeric(tolerance)) {
    stop("'tolerance' must be type numeric.", call. = FALSE)
  }
  
  if (tolerance >= threshold) {
    stop("'tolerance' cannot be >= `threshold`.", call. = FALSE)
  }
  
  if (tolerance < 0 | tolerance > 0.1) {
    stop("'tolerance' must be > 0 and <= 0.1.", call. = FALSE)
  }
  
  if("type" %in% colnames(existing))
  {
    message("`existing` has a variable named `type`. This will cause issues with plotting. Consider removing.")
  }
  

  #--- determine number of bands in mraster ---#

  nb <- terra::nlyr(mraster)

  #--- extract covariates data from mraster ---#

  vals <- terra::as.data.frame(mraster, xy = TRUE, row.names = FALSE) %>%
    dplyr::rename(
      X = x,
      Y = y
    )

  #--- Remove NA / NaN / Inf values ---#

  vals <- vals %>%
    stats::na.omit()

  #--- Generate quantile matrix or use supplied one ---#
  
  if(is.null(matrices)){
    
    #--- if null generate a new quantile matrix ---#
    
    mats <- calculate_pop(mraster = mraster, nQuant = nQuant)
    
  } else {
    
    #--- if quantile matrix is provided ---#

    if (any(!c("matQ", "matCov") %in% names(matrices))) {
      stop("'matrices' must be the output from 'calculate_pop()'.", call. = FALSE)
    }
    
    if(nrow(matrices$matCov) != nQuant){
      stop("Number of quantiles in provided 'matrices' does not match 'nQuant'.", call. = FALSE)
    }
    
    if(!all(names(matrices$values) == names(mraster))) {
      stop("'mraster' used to generate 'matrices' must be identical.", call. = FALSE)
    }
    
    mats <- matrices
    
  }

  #--- Change 0's to very small number to avoid division issues ---#

  mats$matCov[which(mats$matCov == 0)] <- 0.0000001

  #--- Create density matrix from covariates and length of mraster ---#

  matCovDens <- mats$matCov / nrow(vals)

  #--- Remove quantiles that do not cover at least 1% area in each covariate ---#

  matCovDens[which(matCovDens <= 0.01)] <- NA

  ### --- Prepare existing sample data ---###
  if(!inherits(existing, "sf")){
    
    #--- determine crs of input sraster ---#
    
    crs <- terra::crs(mraster, proj = TRUE)
    
    if (any(!c("X", "Y") %in% colnames(existing))) {
      
      #--- if coordinate column names are lowercase change them to uppercase to match requirements ---#
      
      if (any(c("x", "y") %in% colnames(existing))) {
        existing <- existing %>%
          dplyr::rename(
            X = x,
            Y = y
          )
        
        message("'existing' column coordinate names are lowercase - converting to uppercase.")
      } else {
        
        #--- if no x/y columns are present stop ---#
        
        stop("'existing' must have columns named 'X' and 'Y'.", call. = FALSE)
      }
    }
    
    existing <- existing %>%
      sf::st_as_sf(., coords = c("X", "Y"), crs = crs)
    
  } else {
    
    #--- determine crs of input sraster ---#
    
    crs <- sf::st_crs(existing)
    
  }

  #--- extract covariates at existing sample locations ---#

  samples <- extract_metrics(mraster, existing, data.frame = TRUE)
  
  #--- remove already existing samples from vals to not repeat sample ---#
  
  vals <- vals %>%
    dplyr::anti_join(samples, by = c("X", "Y"))
  
  #--- Rearrange columns ---#
  
  #--- do other attributes exist in `existing` - if yes, save them for later ---#
  
  if(length(names(samples))-2 != length(names(mraster))){
    
    extraCols <- samples %>%
      dplyr::select(!names(mraster))
    
  }
  
  #--- Assign attribute to differentiate between original samples and those added during HELS algorithm ---#
  
  samples$type <- "existing"
  
  #--- subset columns for sampling ---#
  
  samples <- samples %>%
    dplyr::select(X, Y, type, names(mraster))
  
  #--- check if samples fall in areas where stratum values are NA ---#
  
  if(any(!complete.cases(samples))){
    
    samples_NA <- samples %>%
      dplyr::filter(!complete.cases(.)) %>%
      dplyr::mutate(type = "existing")

    samples <- samples %>%
      stats::na.omit()
    
  }

  #--- Create data hypercube of existing samples to compare with mraster data ---#

  matCovSamp <- mat_covNB(vals = samples[4:ncol(samples)], nQuant = nQuant, nb = nb, matQ = mats$matQ)

  #--- Change 0's to very small number to avoid division issues ---#

  matCovSamp[which(matCovSamp == 0)] <- 0.0000001

  #--- Create density matrix from covariates and length of mraster ---#

  matCovSampDens <- matCovSamp / nrow(samples)

  ### --- Selection of new samples based on density ---###

  #--- Ratio and ordering of data density and covariate density ---#

  ratioExisting <- ratio <- matCovSampDens / matCovDens

  #--- order the densities based on representation ---#

  #--- low to high ---#

  ratOrderUnder <- order(ratio)

  #--- Outline quantiles that are underrepresented (< 1) in the sample ---#

  underRep <- which(ratio < (1 - tolerance), arr.ind = TRUE)
  underRep <- cbind(underRep, which(ratio < (1 - tolerance)))

  #--- perform sampling ---#

  if (!is.null(nSamp)) {

    out <- ahels_nSamp(nSamp = nSamp,
                       nQuant = nQuant,
                       tolerance = tolerance,
                       nb = nb,
                       underRep = underRep,
                       ratio = ratio,
                       ratOrderUnder = ratOrderUnder,
                       matCovDens = matCovDens,
                       matCovSampDens = matCovSampDens,
                       samples = samples,
                       mats = mats,
                       vals = vals)
    
  } else {
    
    out <- ahels_threshold(threshold = threshold,
                           tolerance = tolerance,
                           ratio = ratio,
                           nQuant = nQuant,
                           nb = nb,
                           underRep = underRep,
                           ratOrderUnder = ratOrderUnder,
                           matCovDens = matCovDens,
                           matCovSampDens = matCovSampDens,
                           samples = samples,
                           mats = mats,
                           vals = vals)

  }

  message(paste0("A total of ", out$sTot, " new samples added."))

  #--- replace existing samples (if they exist) that had NA values for metrics ---#
  
  if(exists("samples_NA")){
    
    if(!is.null(extraCols)){
      
      samples <- out$samples %>%
        dplyr::bind_rows(., samples_NA) %>%
        dplyr::left_join(., extraCols,  by = c("X","Y")) %>%
        sf::st_as_sf(., coords = c("X", "Y"), crs = crs)
      
    } else {
    
      #--- convert coordinates to a spatial points object ---#
      samples <- out$samples %>%
        dplyr::bind_rows(., samples_NA) %>%
        sf::st_as_sf(., coords = c("X", "Y"), crs = crs)
    
    }
    
  } else {
    
    if(!is.null(extraCols)){
      
      samples <- out$samples %>%
        dplyr::left_join(., extraCols,  by = c("X","Y")) %>%
        sf::st_as_sf(., coords = c("X", "Y"), crs = crs)
      
    } else {
      
      #--- convert coordinates to a spatial points object ---#
      samples <- out$samples %>%
        sf::st_as_sf(., coords = c("X", "Y"), crs = crs)
      
    }
    
  }

  #--- assign sraster crs to spatial points object ---#
  sf::st_crs(samples) <- crs

  if (isTRUE(plot)) {
    terra::plot(mraster[[1]])
    suppressWarnings(terra::plot(samples, add = T, col = "black", pch = ifelse(samples$type == "existing", 1, 3)))
  }

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

  #--- output samples & / or samples and details (ratio matrix) ---#

  if (isFALSE(details)) {
    return(samples)
  } else {
    return(list(samples = samples, details = list(
      existingRatio = ratioExisting,
      sampledRatio = out$ratio,
      diffRatio = out$ratio - ratioExisting
    )))
  }
}

