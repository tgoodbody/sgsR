#' CLHS sampling
#'
#' @description Conditioned Latin Hypercube Sampling using the \code{\link[clhs]{clhs-package}} package functionality
#'
#' @family sample functions
#'
#' @inheritParams analyze_sampOptLHC
#' @inheritParams sample_srs
#' @inheritParams strat_kmeans
#' @inheritParams extract_strata
#'
#' @param ... Additional arguments for clhs sampling. See \code{\link[clhs]{clhs}}.
#'
#' @importFrom magrittr %>%
#' @importFrom methods is
#' @importFrom stats coef complete.cases median quantile sd
#' @importFrom utils setTxtProgressBar
#'
#' @return An sf object with \code{nSamp} stratified samples.
#'
#' @export

sample_clhs <- function(mraster,
                        nSamp,
                        existing = NULL,
                        access = NULL,
                        buff_inner = NULL,
                        buff_outer = NULL,
                        plot = FALSE,
                        details = FALSE,
                        ...) {

  #--- check for required packages ---#
  if (!requireNamespace("clhs", quietly = TRUE)) {
    stop("Package \"clhs\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  #--- Set global vars ---#

  x <- y <- type <- NULL

  #--- Error management ---#

  if (!inherits(mraster, "SpatRaster")) {
    stop("'mraster' must be type SpatRaster", call. = FALSE)
  }

  if (!is.numeric(nSamp)) {
    stop("'nSamp' must be type numeric")
  }

  if (is.na(crs(mraster))) {
    stop("'mraster' does not have a coordinate system")
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical")
  }

  if (!is.logical(details)) {
    stop("'details' must be type logical")
  }


  #--- determine crs of input mraster ---#
  crs <- crs(mraster)

  #--- access buffering if specified ---#

  if (!is.null(access)) {

    #--- error handling in the presence of 'access' ---#
    if (!inherits(access, "sf")) {
      stop("'access' must be an 'sf' object")
    }

    if (!inherits(sf::st_geometry(access), "sfc_MULTILINESTRING")) {
      stop("'access' geometry type must be 'sfc_MULTILINESTRING'")
    }

    #--- buffer roads and mask ---#
    
    access_buff <- mask_access(raster = mraster, access = access, buff_inner = buff_inner, buff_outer = buff_outer)
    
    mraster_access <- access_buff$rast

    #--- extract covariate data from mraster ---#

    vals <- terra::as.data.frame(mraster_access, xy = TRUE, row.names = FALSE) %>%
      dplyr::rename(
        X = x,
        Y = y
      )
  } else {

    #--- extract covariate data from mraster ---#

    vals <- terra::as.data.frame(mraster, xy = TRUE, row.names = FALSE) %>%
      dplyr::rename(
        X = x,
        Y = y
      )
  }

  #--- Remove NA / NaN / Inf values ---#

  vals <- vals %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::mutate(type = "new")

  #--- error handing and existing samples preparation ---#

  if (!is.null(existing)) {
    if (!inherits(existing, "data.frame") && !inherits(existing, "sf")) {
      stop("'existing' must be a data.frame or sf object")
    }

    #--- combined existing samples with vals dataframe ---#

    if (!inherits(existing, "sf")) {
      if (any(!c("X", "Y") %in% colnames(existing))) {

        #--- if coordinate column names are lowercase change them to uppercase to match requirements ---#

        if (any(c("x", "y") %in% colnames(existing))) {
          existing <- existing %>%
            dplyr::rename(
              X = x,
              Y = y
            )

          message("Column coordinates names for 'existing' are lowercase - converting to uppercase")
        } else {

          #--- if no x/y columns are present stop ---#

          stop("'existing' must have columns named 'X' and 'Y'")
        }
      }

      existingSamples <- existing
    } else {
      existingSamples <- extract_metrics(mraster, existing, data.frame = TRUE)
    }

    #--- create dataset with labels for plotting ---#

    existingSamples <- existingSamples %>%
      dplyr::mutate(type = "existing")

    #--- create conjoined existing dataset ---#

    vals <- rbind(existingSamples, vals)
  }

  vals_tp <- vals %>% dplyr::select(-type)

  ##########################
  #### SAMPLING ############
  ##########################

  #--- if existing samples are not provided ---#

  if (is.null(existing)) {

    #--- remove 'type' during sampling ---#

    clhsOut <- clhs::clhs(x = vals_tp, size = nSamp, ...)


    #--- if ... variables are provided the output is sometimes a list object ---#

    if (is.list(clhsOut)) {
      samples <- clhsOut$sampled_data
    } else {

      #--- take samples from original vals dataframe for 'type' attribute ---#

      samples <- vals[clhsOut, ]
    }
  } else {

    #--- same as above but this time including existing samples ---#

    clhsOut <- clhs::clhs(x = vals_tp, size = nSamp, include = 1:nrow(existingSamples), ...)

    if (inherits(clhsOut, "list")) {

      #--- locate row indices for each sample and extract ---#

      index <- clhsOut$index_samples

      samples <- vals[index, ]
    } else {
      samples <- vals[clhsOut, ]
    }
  }

  #--- convert coordinates to a spatial points object ---#
  samples <- samples %>%
    as.data.frame() %>%
    sf::st_as_sf(., coords = c("X", "Y"))

  #--- assign sraster crs to spatial points object ---#
  sf::st_crs(samples) <- crs

  if (isTRUE(plot)) {
    if (!is.null(access)) {

      #--- plot samples as well as aggregated access buffers ---#

      terra::plot(mraster[[1]])
      suppressWarnings(terra::plot(access_buff$buff, add = T, border = c("gray30"), col = "gray10", alpha = 0.1))
      suppressWarnings(terra::plot(samples, add = T, col = ifelse(samples$type == "existing", "Black", "Red")))
    } else {
      terra::plot(mraster[[1]])
      suppressWarnings(terra::plot(samples, add = T, col = "black", pch = ifelse(samples$type == "existing", 1, 3)))
    }
  }

  if (isTRUE(details)) {

    #--- output metrics details along with stratification raster ---#

    output <- list(clhs = clhsOut, samples = samples)

    #--- output samples dataframe ---#
    return(output)
  } else {

    #--- just output raster ---#

    return(samples)
  }
}
