#' Adapted Hypercube Evaluation of a Legacy Sample (ahels)
#'
#' @description Perform the adapted Hypercube Evaluation of a Legacy Sample (ahels) algorithm using
#' existing site data and raster metrics. New samples are allocated based on quantile ratios between
#' the existing sample and covariate dataset.
#'
#' @family sample functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams extract_strata
#'
#' @param nQuant Numeric. Number of quantiles to divide covariates and samples into. Quantiles that do not
#' cover at least 1 percent of the area of interest will be excluded and be returned as \code{NA}.
#' @param nSamp Numeric. Maximum number of new samples to allocate. If provided, the algorithm
#' will default to allocating number of samples provided.
#' @param threshold Numeric. A sample quantile ratio threshold for establishing whether
#' additional samples should be added. \code{default = 0.9}. Values close to 1 can cause the algorithm to
#' continually loop and should be used sparingly.
#' @param matCov List. Covariance matrix generated from \code{calculate_lhsPop(mraster = mraster, PCA = FALSE, nQuant = nQuant)}.
#' Both \code{mraster} & \code{nQuant} inputs must be the same to supply the covariance matrix. Supplying the matrix allows users
#' with very large rasters to pre-process the covariance matrix to avoid longer sampling processing times.
#' @param plot Logial. Plots existing (circles) and new (crosses) samples on the first band of mraster.
#'
#' @references
#' Malone BP, Minansy B, Brungard C. 2019. Some methods to improve the utility of conditioned Latin hypercube sampling. PeerJ 7:e6451 DOI 10.7717/peerj.6451
#'
#' @return Returns sf point object with existing samples and supplemental samples added by the ahels algorithm.
#'
#' @examples
#' #--- Load raster and existing plots---#
#' r <- system.file("extdata", "wall_metrics.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' e <- system.file("extdata", "existing.shp", package = "sgsR")
#' e <- sf::st_read(e)
#'
#' sample_ahels(
#'   mraster = mr[[1:3]],
#'   existing = e,
#'   plot = TRUE
#' )
#'
#' sample_ahels(
#'   mraster = mr[[1:3]],
#'   existing = e,
#'   nQuant = 20,
#'   nSamp = 300,
#'   filename = tempfile(fileext = ".shp")
#' )
#' @note
#'
#' Messages in the algorithm will state that samples have been added to under-represented quantiles. The number between
#' square brackets that follow represent the matrix row and column respectively that can be printed using \code{details = TRUE}.
#' 
#' In some cases, generally when a single metric is used as \code{mraster}, sampling ratios all be >= 1 before the 
#' \code{nSamp} number of samples are allocated. The algorithm will stop in this case.
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
                         matCov = NULL,
                         plot = FALSE,
                         details = FALSE,
                         filename = NULL,
                         overwrite = FALSE) {

  #--- Set global vars ---#

  x <- y <- X <- Y <- n <- type <- geometry <- NULL

  #--- Error handling ---#

  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  }

  if (!inherits(existing, "data.frame") && !inherits(existing, "sf")) {
    stop("'existing' must be a data.frame or sf object", call. = FALSE)
  }

  if (!is.numeric(nQuant)) {
    stop("'nQuant' must be type numeric", call. = FALSE)
  }

  if (!is.numeric(threshold)) {
    stop("'threshold' must be type numeric", call. = FALSE)
  }

  if (threshold < 0 | threshold > 1) {
    stop("'threshold' must be > 0 and < 1", call. = FALSE)
  }

  if (!is.logical(details)) {
    stop("'details' must be type logical", call. = FALSE)
  }

  #--- determine number of bands in mraster ---#

  nb <- terra::nlyr(mraster)

  #--- determine crs of input sraster ---#

  crs <- terra::crs(mraster, proj = TRUE)

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
  
  if(is.null(matCov)){
    
    #--- if null generate a new quantile matrix ---#
    
    mats <- calculate_lhsPop(mraster = mraster, PCA = FALSE, nQuant = nQuant)
    
  } else {
    
    #--- if quantile matrix is provided ---#
    
    if (any(!c("values","pcaLoad","matQ", "matCov") %in% names(matCov))) {
      stop("'matCov' must be the output from `calculate_lhsPop()`.", call. = FALSE)
    }
    
    if(nrow(matCov$matCov) > nQuant){
      stop("Number of quantiles in provided 'matCov' does not match nQuant.", call. = FALSE)
    }
    
    if(any(!names(matCov$values) %in% names(mraster))){
      message("'mraster' used to generate 'matCov' must be identical,", call. = FALSE)
    }
    
    mats <- matCov
    
  }

  #--- Change 0's to very small number to avoid division issues ---#

  mats$matCov[which(mats$matCov == 0)] <- 0.0000001

  #--- Create density matrix from covariates and length of mraster ---#

  matCovDens <- mats$matCov / nrow(vals)

  #--- Remove quantiles that do not cover at least 1% area in each covariate ---#

  matCovDens[which(matCovDens <= 0.01)] <- NA

  ### --- Prepare existing sample data ---###

  #--- select geometry attribute ---#

  existing <- existing %>%
    dplyr::select(geometry)

  #--- extract covariates at existing sample locations ---#

  samples <- extract_metrics(mraster, existing, data.frame = TRUE)

  #--- remove already existing samples from vals to not repeat sample ---#

  vals <- vals %>%
    dplyr::anti_join(samples, by = c("X", "Y"))

  #--- Assign attribute to differentiate between original samples and those added during HELS algorithm ---#

  samples$type <- "existing"
  samples$n <- seq(1:nrow(samples))

  #--- Rearrange columns ---#

  samples <- samples %>%
    dplyr::select(X, Y, n, type, tidyselect::everything())

  #--- Create data hypercube of existing samples to compare with mraster data ---#

  matCovSamp <- mat_covNB(vals = samples[5:ncol(samples)], nQuant = nQuant, nb = nb, matQ = mats$matQ)

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

  underRep <- which(ratio < 1, arr.ind = TRUE)
  underRep <- cbind(underRep, which(ratio < 1))

  #--- begin sampling from highest discrepancy to lowest ---#

  newSamp <- nSamp
  sTot <- 0

  #--- begin while loop to sample ---#


  if (!is.null(nSamp)) {

    #--- ensure nSamp is numeric ---#

    if (!is.numeric(nSamp)) {
      stop("'nSamp' must be type numeric", call. = FALSE)
    }

    message(glue::glue("nSamp of {nSamp} has been provided. Samples will be added until this number is reached or until sampling ratios are all >= 1"))

    while (newSamp != 0) {

      #--- determine the greatest discrepancy between sample and covariate data ---#

      repRankUnder <- which(underRep[, 3] == ratOrderUnder[1])

      #--- determine row and column of most under represented quantile ---#

      repRow <- underRep[repRankUnder, 1]
      repCol <- underRep[repRankUnder, 2]
      
      #--- if all sampling ratios in matCovSampDens are >= 1 stop adding samples ---#
      
      if(length(repRankUnder) == 0){
        
        message(glue::glue("Sampling ratios are all >= 1. A total of {sTot} samples were added."))
        
        break
      }

      #--- determine number of existing samples in selected quantile ---#
      sampExist <- floor(matCovSamp[repRow, repCol])

      #--- determine max number of samples based on covariate density ---#

      sampOptim <- ceiling(nrow(samples) * matCovDens[repRow, repCol])

      #--- number of samples needed ---#

      sampNeed <- sampOptim - sampExist
      
      #--- we have a limited number of samples so we need to be sure not to over allocate ---#
      if (newSamp <= sampNeed) sampNeed <- newSamp
      
      #--- selecting covariates based on quantile chosen ---#

      covLower <- mats$matQ[repRow, repCol]

      covUpper <- mats$matQ[repRow + 1, repCol]

      #--- subset covariate dataset for potential new samples ---#

      valsSub <- vals[vals[, (2 + repCol)] >= covLower & vals[, (2 + repCol)] <= covUpper, ]

      #--- randomly sample within valsSub and extract randomly sampled cells ---#

      addSamp <- sample(nrow(valsSub), sampNeed)

      valsSubSamp <- valsSub[addSamp, ]

      valsSubSamp$type <- "new"
      valsSubSamp$n <- row.names(valsSubSamp)

      #--- remove samples from pool to ensure same cells are not sampled again ---#

      vals <- vals[-addSamp, ]

      #--- add new samples to existing sample dataframe ---#

      samples <- rbind(samples, valsSubSamp)

      #--- update loop parameters ---#

      message("Under-represented Quantile ", paste0("[",repRow, ",",repCol, "]"), " - A total of ", sampNeed, " samples have been allocated.")

      #--- update total allocated samples ---#

      sTot <- sTot + sampNeed

      #--- update available sample number ---#

      newSamp <- newSamp - sampNeed

      #--- recompute ratio's in the presence of newly added samples ---#

      matCovSamp <- mat_covNB(vals = samples[5:ncol(samples)], nQuant = nQuant, nb = nb, matQ = mats$matQ)

      #--- Change 0's to very small number to avoid division issues ---#

      matCovSamp[which(matCovSamp == 0)] <- 0.0000001

      #--- Create density matrix from covariates and length of mraster ---#

      matCovSampDens <- matCovSamp / nrow(samples)

      ### --- Selection of new samples based on density ---###

      #--- Ratio and ordering of data density and covariate density ---#

      ratio <- matCovSampDens / matCovDens

      #--- order the densities based on representation ---#

      #--- low to high ---#

      ratOrderUnder <- order(ratio)

      #--- Outline quantiles that are underrepresented (< 1) in the sample ---#

      underRep <- which(ratio < 1, arr.ind = TRUE)
      underRep <- cbind(underRep, which(ratio < 1))
    }
  } else {
    message(glue::glue("Threshold of {threshold} has been provided. Samples will be added until sampling ratios are >= {threshold}."))

    ### --- If 'nSamp' is not provided a threshold is used ---###

    while (isTRUE(any(ratio < threshold))) {

      #--- determine the greatest discrepancy between sample and covariate data ---#

      repRankUnder <- which(underRep[, 3] == ratOrderUnder[1])

      #--- determine row and column of most under represented quantile ---#

      repRow <- underRep[repRankUnder, 1]
      repCol <- underRep[repRankUnder, 2]

      #--- determine number of existing samples in selected quantile ---#

      sampExist <- floor(nrow(samples) * matCovSampDens[repRow, repCol])

      #--- determine max number of samples based on covariate density ---#

      sampOptim <- ceiling(nrow(samples) * matCovDens[repRow, repCol])

      #--- number of samples needed ---#

      sampNeed <- sampOptim - sampExist

      #--- selecting covariates based on quantile chosen ---#

      covLower <- mats$matQ[repRow, repCol]

      covUpper <- mats$matQ[repRow + 1, repCol]

      #--- subset covariate dataset for potential new samples ---#

      valsSub <- vals[vals[, (2 + repCol)] >= covLower & vals[, (2 + repCol)] <= covUpper, ]

      #--- randomly sample within valsSub and extract randomly sampled cells ---#

      addSamp <- sample(nrow(valsSub), sampNeed)

      valsSubSamp <- valsSub[addSamp, ]

      valsSubSamp$type <- "new"

      valsSubSamp$n <- row.names(valsSubSamp)

      #--- remove samples from pool to ensure same cells are not sampled again ---#

      vals <- vals[-addSamp, ]

      #--- add new samples to existing sample dataframe ---#

      samples <- rbind(samples, valsSubSamp)

      #--- update loop parameters ---#

      message("Under-represented Quantile ", paste0("[",repRow, ",",repCol, "]"), " - A total of ", sampNeed, " samples have been allocated.")

      #--- update total allocated samples ---#

      sTot <- sTot + sampNeed

      #--- recompute ratio's in the presence of newly added samples ---#

      matCovSamp <- mat_covNB(vals = samples[5:ncol(samples)], nQuant = nQuant, nb = nb, matQ = mats$matQ)

      #--- Change 0's to very small number to avoid division issues ---#

      matCovSamp[which(matCovSamp == 0)] <- 0.0000001

      #--- Create density matrix from covariates and length of mraster ---#

      matCovSampDens <- matCovSamp / nrow(samples)

      ### --- Selection of new samples based on density ---###

      #--- Ratio and ordering of data density and covariate density ---#

      ratio <- matCovSampDens / matCovDens

      #--- order the densities based on representation ---#

      #--- low to high ---#

      ratOrderUnder <- order(ratio)

      #--- high to low ---#

      ratOrderOver <- rev(ratOrderUnder)

      #--- Outline quantiles that are underrepresented (< 1) in the sample ---#

      underRep <- which(ratio < 1, arr.ind = TRUE)
      underRep <- cbind(underRep, which(ratio < 1))
    }
  }

  message(glue::glue("A total of {sTot} new samples added"))

  #--- convert coordinates to a spatial points object ---#
  samples <- samples %>%
    as.data.frame() %>%
    dplyr::select(-n) %>%
    sf::st_as_sf(., coords = c("X", "Y"))

  #--- assign sraster crs to spatial points object ---#
  sf::st_crs(samples) <- crs

  if (isTRUE(plot)) {
    terra::plot(mraster[[1]])
    suppressWarnings(terra::plot(samples, add = T, col = "black", pch = ifelse(samples$type == "existing", 1, 3)))
  }

  if (!is.null(filename)) {
    if (!is.logical(overwrite)) {
      stop("'overwrite' must be either TRUE or FALSE")
    }

    if (file.exists(filename) & isFALSE(overwrite)) {
      stop(glue::glue("{filename} already exists and overwrite = FALSE"))
    }

    sf::st_write(samples, filename, delete_layer = overwrite)
  }

  #--- output samples & / or samples and details (ratio matrix) ---#

  if (isFALSE(details)) {
    return(samples)
  } else {
    return(list(samples = samples, details = list(
      existingRatio = ratioExisting,
      sampledRatio = ratio,
      diffRatio = ratio - ratioExisting
    )))
  }
}

###--- nSamp ahels algorithm ---###
ahels_nSamp <- function(newSamp,
                        sTot,
                        underRep,
                        ratOrderUnder,
                        mats,
                        vals,
                        matCovSamp,
                        matCovDens,
                        samples,
                        ){
  
  while (newSamp != 0) {
    
    #--- determine the greatest discrepancy between sample and covariate data ---#
    
    repRankUnder <- which(underRep[, 3] == ratOrderUnder[1])
    
    #--- determine row and column of most under represented quantile ---#
    
    repRow <- underRep[repRankUnder, 1]
    repCol <- underRep[repRankUnder, 2]
    
    #--- if all sampling ratios in matCovSampDens are >= 1 stop adding samples ---#
    
    if(length(repRankUnder) == 0){
      
      message(glue::glue("Sampling ratios are all >= 1. A total of {sTot} samples were added."))
      
      break
    }
    
    #--- determine number of existing samples in selected quantile ---#
    sampExist <- floor(matCovSamp[repRow, repCol])
    
    #--- determine max number of samples based on covariate density ---#
    
    sampOptim <- ceiling(nrow(samples) * matCovDens[repRow, repCol])
    
    #--- number of samples needed ---#
    
    sampNeed <- sampOptim - sampExist
    
    #--- we have a limited number of samples so we need to be sure not to over allocate ---#
    if (newSamp <= sampNeed) sampNeed <- newSamp
    
    #--- selecting covariates based on quantile chosen ---#
    
    covLower <- mats$matQ[repRow, repCol]
    
    covUpper <- mats$matQ[repRow + 1, repCol]
    
    #--- subset covariate dataset for potential new samples ---#
    
    valsSub <- vals[vals[, (2 + repCol)] >= covLower & vals[, (2 + repCol)] <= covUpper, ]
    
    #--- randomly sample within valsSub and extract randomly sampled cells ---#
    
    addSamp <- sample(nrow(valsSub), sampNeed)
    
    valsSubSamp <- valsSub[addSamp, ]
    
    valsSubSamp$type <- "new"
    valsSubSamp$n <- row.names(valsSubSamp)
    
    #--- remove samples from pool to ensure same cells are not sampled again ---#
    
    vals <- vals[-addSamp, ]
    
    #--- add new samples to existing sample dataframe ---#
    
    samples <- rbind(samples, valsSubSamp)
    
    #--- update loop parameters ---#
    
    message("Under-represented Quantile ", paste0("[",repRow, ",",repCol, "]"), " - A total of ", sampNeed, " samples have been allocated.")
    
    #--- update total allocated samples ---#
    
    sTot <- sTot + sampNeed
    
    #--- update available sample number ---#
    
    newSamp <- newSamp - sampNeed
    
    #--- recompute ratio's in the presence of newly added samples ---#
    
    matCovSamp <- mat_covNB(vals = samples[5:ncol(samples)], nQuant = nQuant, nb = nb, matQ = mats$matQ)
    
    #--- Change 0's to very small number to avoid division issues ---#
    
    matCovSamp[which(matCovSamp == 0)] <- 0.0000001
    
    #--- Create density matrix from covariates and length of mraster ---#
    
    matCovSampDens <- matCovSamp / nrow(samples)
    
    ### --- Selection of new samples based on density ---###
    
    #--- Ratio and ordering of data density and covariate density ---#
    
    ratio <- matCovSampDens / matCovDens
    
    #--- order the densities based on representation ---#
    
    #--- low to high ---#
    
    ratOrderUnder <- order(ratio)
    
    #--- Outline quantiles that are underrepresented (< 1) in the sample ---#
    
    underRep <- which(ratio < 1, arr.ind = TRUE)
    underRep <- cbind(underRep, which(ratio < 1))
  }
  
  return(samples)
}
