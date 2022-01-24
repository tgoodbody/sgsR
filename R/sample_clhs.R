#' CLHS sampling
#'
#' @description Conditioned Latin Hypercube Sampling using \code{\link[clhs]{clhs}} functionality.
#'
#' @family sample functions
#'
#' @inheritParams calculate_lhsOpt
#' @inheritParams sample_srs
#' @inheritParams strat_kmeans
#' @inheritParams extract_strata
#'
#' @param iter Numeric. Value giving the number of iterations within the Metropolis-Hastings process.
#' @param cost Numeric/Character. Index or name of covariate within \code{mraster} to be used to constrain cLHS sampling.
#' If default - \code{NULL} then a cost constraint is not used.
#' @param ... Additional arguments for clhs sampling. See \code{\link[clhs]{clhs}}.
#'
#' @importFrom stats coef complete.cases median quantile sd
#' @importFrom utils setTxtProgressBar
#'
#' @return An sf object with \code{nSamp} stratified samples.
#'
#' @examples
#' #--- Load raster and existing plots---#
#' r <- system.file("extdata", "wall_metrics.tif", package = "sgsR")
#' mr <- terra::rast(r)
#'
#' e <- system.file("extdata", "existing.shp", package = "sgsR")
#' e <- sf::st_read(e)
#'
#' a <- system.file("extdata", "roads.shp", package = "sgsR")
#' ac <- sf::st_read(a)
#'
#' sample_clhs(
#'   mraster = mr,
#'   nSamp = 200,
#'   plot = TRUE,
#'   iter = 100
#' )
#'
#' sample_clhs(
#'   mraster = mr,
#'   nSamp = 400,
#'   existing = e,
#'   iter = 250,
#'   details = TRUE
#' )
#'
#' sample_clhs(
#'   mraster = mr,
#'   nSamp = 200,
#'   iter = 200,
#'   existing = e,
#'   access = ac,
#'   buff_inner = 100,
#'   buff_outer = 300,
#'   plot = TRUE
#' )
#'
#' #--- cost constrained examples ---#
#' #--- calculate distance to access layer for each pixel in mr ---#
#' mr.c <- calculate_distance(
#'   raster = mr,
#'   access = ac
#' )
#'
#' sample_clhs(
#'   mraster = mr.c,
#'   nSamp = 250,
#'   iter = 200,
#'   cost = "dist2access",
#'   plot = TRUE
#' )
#'
#' sample_clhs(
#'   mraster = mr.c,
#'   nSamp = 250,
#'   existing = e,
#'   iter = 200,
#'   cost = "dist2access",
#'   plot = TRUE
#' )
#' @references
#' Minasny, B. and McBratney, A.B. 2006. A conditioned Latin hypercube method
#' for sampling in the presence of ancillary information. Computers and
#' Geosciences, 32:1378-1388.
#'
#' Minasny, B. and A. B. McBratney, A.B.. 2010. Conditioned Latin Hypercube
#' Sampling for Calibrating Soil Sensor Data to Soil Properties. In: Proximal
#' Soil Sensing, Progress in Soil Science, pages 111-119.
#'
#' Roudier, P., Beaudette, D.E. and Hewitt, A.E. 2012. A conditioned Latin
#' hypercube sampling algorithm incorporating operational constraints. In:
#' Digital Soil Assessments and Beyond. Proceedings of the 5th Global Workshop
#' on Digital Soil Mapping, Sydney, Australia.
#'
#' @author Tristan R.H. Goodbody
#'
#' @export

sample_clhs <- function(mraster,
                        nSamp,
                        iter = 10000,
                        cost = NULL,
                        existing = NULL,
                        access = NULL,
                        buff_inner = NULL,
                        buff_outer = NULL,
                        plot = FALSE,
                        details = FALSE,
                        filename = NULL,
                        overwrite = FALSE,
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

  if (!is.numeric(iter)) {
    stop("'iter' must be type numeric")
  }

  if (is.na(terra::crs(mraster, proj = TRUE))) {
    stop("'mraster' does not have a coordinate system")
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical")
  }

  if (!is.logical(details)) {
    stop("'details' must be type logical")
  }


  #--- determine crs of input mraster ---#
  crs <- terra::crs(mraster, proj = TRUE)

  #--- access buffering if specified ---#

  if (!is.null(access)) {

    #--- error handling in the presence of 'access' ---#
    if (!inherits(access, "sf")) {
      stop("'access' must be an 'sf' object")
    }

    if (!inherits(sf::st_geometry(access), "sfc_MULTILINESTRING") && !inherits(sf::st_geometry(access), "sfc_LINESTRING")) {
      stop("'access' geometry type must be 'LINESTRING' or 'MULTILINESTRING'")
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

  #--- incorporate cost constraint ---#

  if (!is.null(cost)) {
    if (is.numeric(cost)) {
      if ((cost + 2) > (ncol(vals)) | cost < 0) {
        stop("'cost' index doest not exist within 'mraster'")
      }

      #--- need to add 2 because X and Y are added to vals ---#

      cost <- cost + 2
    } else if (is.character(cost)) {
      cost <- which(names(vals) == cost)
    } else {
      stop("'cost' must be either a numeric index or name of covariate within 'mraster'")
    }
  }

  #--- Remove NA / NaN / Inf values ---#

  vals <- vals %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::mutate(type = "new")

  namesvals <- names(vals)

  #--- error handing and existing samples preparation ---#

  if (!is.null(existing)) {
    if (!inherits(existing, "data.frame") && !inherits(existing, "sf")) {
      stop("'existing' must be a data.frame or sf object")
    }

    #--- check that nSamp is > than existing ---#

    if (nrow(existing) > nSamp) {
      stop("nSamp must be > than number of existing samples")
    }

    #--- combine existing samples with vals dataframe ---#

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
      dplyr::mutate(type = "existing") %>%
      dplyr::select(namesvals)

    #--- create conjoined existing dataset ---#

    vals <- rbind(existingSamples, vals)
  }

  #--- remove 'type' during sampling ---#
  
  vals_tp <- vals %>% dplyr::select(-type)

  ##########################
  #### SAMPLING ############
  ##########################

  #--- if existing samples are not provided ---#

  if (is.null(existing)) {

    clhsOut <- clhs::clhs(x = vals_tp, size = nSamp, iter = iter, cost = cost, ...)

    #--- if ... variables are provided the output is sometimes a list object ---#

    if (is.list(clhsOut)) {
      samples <- clhsOut$sampled_data
    } else {

      #--- take samples from original vals dataframe for 'type' attribute ---#

      samples <- vals[clhsOut, ]
    }
  } else {

    #--- same as above but this time including existing samples ---#

    clhsOut <- clhs::clhs(x = vals_tp, size = nSamp, iter = iter, cost = cost, include = 1:nrow(existingSamples), ...)

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
      suppressWarnings(terra::plot(samples, add = T, col = "black", pch = ifelse(samples$type == "existing", 1, 3)))
    } else {
      terra::plot(mraster[[1]])
      suppressWarnings(terra::plot(samples, add = T, col = "black", pch = ifelse(samples$type == "existing", 1, 3)))
    }
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
