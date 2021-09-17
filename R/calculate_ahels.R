#' Adapted Hypercube Evaluation of a Legacy Sample (ahels)
#'
#' @description Perform the adapted Hypercube Evaluation of a Legacy Sample (ahels) algorithm using
#' existing site data and raster metrics. New samples are allocated based on quantile ratios between
#' the existing sampleand covariate dataset.
#'
#' @family sample functions
#'
#' @inheritParams strat_kmeans
#' @inheritParams extract_strata
#'
#' @param nSamp Numeric. Maximum number of new samples to allocate. If provided, the algorithm
#' will default to allocating number of samples provided.
#' @param threshold Numeric. A sample quantile ratio threshold for establishing whether
#' additional samples should be added. \code{default = 0.9}. Values close to 1 can cause the algorithm to
#' continually loop and should be used sparingly.
#' @param plot Logial. Plots existing (circles) and new (crosses) samples on the first band of mraster.
#' @param nQuant Numeric. Number of quantiles to divide covariates and samples into.
#'
#' @references
#' Malone BP, Minansy B, Brungard C. 2019. Some methods to improve the utility of conditioned Latin hypercube sampling. PeerJ 7:e6451 DOI 10.7717/peerj.6451
#'
#' @return Returns sf point object with existing samples and supplemental samples added by the ahels algorithm.
#'
#' @examples
#' #--- Load raster and existing plots---#
#' r <- system.file("extdata", "wall_metrics_small.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' e <- system.file("extdata", "existing.shp", package = "sgsR")
#' e <- sf::st_read(e)
#' 
#' sample_ahels(mraster = mr[[1:3]], 
#'              existing = e, 
#'              plot = TRUE)
#' 
#' sample_ahels(mraster = mr[[1:3]], 
#'              existing = e, 
#'              nQuant = 20, 
#'              nSamp = 300,
#'              filename = tempfile(fileext = ".shp"))
#' @note
#' Special thanks to Brendan Malone for the original implementation of this algorithm.
#'
#' @author Tristan R.H. Goodbody
#'
#' @export

sample_ahels <- function(mraster,
                         existing,
                         nQuant = 10,
                         nSamp = NULL,
                         threshold = 0.9,
                         plot = FALSE,
                         filename = NULL,
                         overwrite = FALSE) {

  #--- Set global vars ---#

  x <- y <- X <- Y <- n <- type <- geometry <- NULL

  #--- Error handling ---#

  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  }

  if (!inherits(existing, "data.frame") && !inherits(existing, "sf")) {
    stop("'existing' must be a data.frame or sf object")
  }

  if (!is.numeric(nQuant)) {
    stop("'nQuant' must be type numeric")
  }

  if (!is.numeric(threshold)) {
    stop("'threshold' must be type numeric")
  }

  if (threshold < 0 | threshold > 1) {
    stop("'threshold' must be > 0 and < 1")
  }

  #--- determine number of bands in 'mraster' ---#

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
    dplyr::filter(complete.cases(.))

  #--- Generate quantile matrix ---#

  mats <- calculate_lhsPop(mraster = mraster, PCA = FALSE, nQuant = nQuant)

  #--- Change 0's to very small number to avoid division issues ---#

  mats$matCov[which(mats$matCov == 0)] <- 0.0000001

  #--- Create density matrix from covariates and length of mraster ---#

  matCovDens <- mats$matCov / nrow(vals)

  #--- Remove quantiles that do not cover at least 1% area in eac covariate ---#

  matCovDens[which(matCovDens <= 0.01)] <- NA

  ### --- Prepare existing sample data ---###

  #--- remove any attributes that are not geometry ---#

  existing <- existing %>%
    dplyr::select(geometry)

  #--- extract covariates at existing sample locations ---#

  samples <- extract_metrics(mraster, existing, data.frame = TRUE)

  #--- remove already existing samples from vals to no repeat sample ---#

  vals <- vals %>%
    dplyr::anti_join(samples, by = c("X", "Y"))

  #--- Assign code to differentiate between original samples and those added during HELS algorithm ---#

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

  ratio <- matCovSampDens / matCovDens

  #--- order the densities based on representation ---#

  #--- low to high ---#

  ratOrderUnder <- order(ratio)

  #--- high to low ---#

  ratOrderOver <- rev(ratOrderUnder)

  #--- Outline quantiles that are underrepresented (< 1) in the sample ---#

  underRep <- which(ratio < 1, arr.ind = TRUE)
  underRep <- cbind(underRep, which(ratio < 1))

  #--- begin sampling from highest discrepancy to lowest ---#

  newSamp <- nSamp
  nLoop <- 1
  sTot <- 0

  #--- begin while loop to sample ---#


  if (!is.null(nSamp)) {

    #--- ensure nSamp is numeric ---#

    if (!is.numeric(nSamp)) {
      stop("'nSamp' must be type numeric")
    }

    message(glue::glue('nSamp of {nSamp} has been provided. Samples will be added until this number is reached'))

    while (newSamp != 0) {

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

      message("Underrepresented Quantile ", nLoop, " - A total of ", sampNeed, " samples have been allocated.")

      #--- update loop number ---#

      nLoop <- nLoop + 1

      #--- update total allocated samples ---#

      sTot <- sTot + sampNeed

      #--- update available sample number ---#

      newSamp <- newSamp - sampNeed

      #---
      #--- recompute ratio's in the presence of newly added samples ---#
      #---

      matCovSamp <- mat_covNB(vals = samples[5:ncol(samples)], nQuant = nQuant, nb = nb, matQ = mats$matQ)

      #--- Change 0's to very small number to avoid division issues ---#

      matCovSamp[which(matCovSamp == 0)] <- 0.0000001

      #--- Create density matrix from covariates and length of mraster ---#

      matCovSampDens <- matCovSamp / nrow(samples)

      ### --- Selection of new samples based on density ---###

      #--- Ratio and ordering of data density and covariate density ---#

      ratio <- matCovSampDens / matCovDens

      # print(ratio)

      #--- order the densities based on representation ---#

      #--- low to high ---#

      ratOrderUnder <- order(ratio)

      #--- high to low ---#

      ratOrderOver <- rev(ratOrderUnder)

      #--- Outline quantiles that are underrepresented (< 1) in the sample ---#

      underRep <- which(ratio < 1, arr.ind = TRUE)
      underRep <- cbind(underRep, which(ratio < 1))
    }
  } else {
    message(glue::glue('threshold of {threshold} has been provided. Samples will be added until quantile ratio is reached'))

    #---
    #--- If 'nSamp' is not provided a threshold is used ---#
    #---

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

      message("Underrepresented Quantile ", nLoop, " - A total of ", sampNeed, " samples have been allocated.")

      #--- update loop number ---#

      nLoop <- nLoop + 1

      #--- update total allocated samples ---#

      sTot <- sTot + sampNeed

      #--- update available sample number ---#

      # newSamp <- newSamp - sampNeed

      #---
      #--- recompute ratio's in the presence of newly added samples ---#
      #---

      matCovSamp <- mat_covNB(vals = samples[5:ncol(samples)], nQuant = nQuant, nb = nb, matQ = mats$matQ)

      #--- Change 0's to very small number to avoid division issues ---#

      matCovSamp[which(matCovSamp == 0)] <- 0.0000001

      #--- Create density matrix from covariates and length of mraster ---#

      matCovSampDens <- matCovSamp / nrow(samples)

      ### --- Selection of new samples based on density ---###

      #--- Ratio and ordering of data density and covariate density ---#

      ratio <- matCovSampDens / matCovDens

      # print(ratio)

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

  message(glue::glue('A total of {sTot} new samples added'))

  print(ratio)

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
      stop(glue::glue('{filename} already exists and overwrite = FALSE'))
    }

    sf::st_write(samples, filename, delete_layer = overwrite)
  }

  return(samples)
}
