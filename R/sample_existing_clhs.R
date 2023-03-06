#' Sub-sample using the conditional Latin hypercube sampling (CLHS)
#'
#' This function takes an existing \code{sf} object and returns a random sub-sample of size \code{nSamp}.
#'
#' @inheritParams sample_clhs
#' @inheritParams extract_strata
#' @inheritParams sample_existing_balanced
#' 
#' @return A subsampled SpatialPointsDataFrame object.
#' 
#' @keywords internal
sample_existing_clhs <- function(existing,
                                 nSamp,
                                 raster = NULL,
                                 cost = NULL,
                                 iter = 10000,
                                 details = FALSE,
                                 filename = NULL,
                                 overwrite = FALSE,
                                 ...) {
  x <- y <- pecdfcat <- pecdf <- NULL
  #--- incorporate cost constraint ---#

  if (!is.null(raster)) {
    if (length(names(raster)) <= 1) {
      stop("At least 2 raster attributes are required to generate a matrix for sub-sampling.", call. = FALSE)
    }
  }

  if (!is.numeric(iter)) {
    stop("'iter' must be type numeric.", call. = FALSE)
  }

  if (iter <= 0) {
    stop("'iter' must be  >= 0.", call. = FALSE)
  }

  #--- get existing values ---#
  vals <- coords_existing(existing)

  if (!is.null(cost)) {
    if (!is.character(cost) & !is.numeric(cost)) {
      stop("'cost' must be either type numeric or character.", call. = FALSE)
    }

    if (is.null(raster)) {
      costLoc <- existing

      v <- "existing"
    } else {
      costLoc <- raster

      v <- "raster"
    }

    if (is.numeric(cost)) {
      if (cost > length(names(costLoc)) | cost < 0) {
        stop(paste0("'cost' index does not exist within '", v, "'."), call. = FALSE)
      }
    } else {
      if (length(which(names(costLoc) == cost)) == 0) {
        stop(paste0("No layer named '", cost, "' exists in '", v, "'."), call. = FALSE)
      } else {
        cost <- which(names(costLoc) == cost)
      }
    }

    message(paste0("Using `", names(costLoc)[cost], "` as sampling constraint."))
  } else {
    costLoc <- NULL
  }

  #--- determine if sampling will be performed using metric attributes in 'existing' only or with the addition of 'raster' metrics ---#

  if (is.null(raster)) {
    message("Sub-sampling based on ALL 'existing' metric distributions. Ensure only attributes of interest are included.")

    all <- existing %>%
      sf::st_drop_geometry(.)

    #--- test if existing has attributes other than geometry ---#
    if (ncol(all) <= 1) {
      # if(names(existing) == "geometry"){
      stop("At least 2 attributes are required to generate a matrix for sub-sampling.", call. = FALSE)
      # }
    }

    #--- sampling ---#

    sidx <- (nrow(all) - nrow(existing) + 1):nrow(all)

    if (isTRUE(details)) {
      #--- output clhs information to be supplied in 'details' list output ---#

      clhsOut <- clhs::clhs(x = all, size = nSamp, iter = iter, cost = cost, can.include = sidx, simple = FALSE, ...)

      outIdx <- clhsOut$index_samples
      #--- extract sampled rows from existing ---#

      samples <- existing[outIdx, ]
    } else {
      outIdx <- clhs::clhs(x = all, size = nSamp, iter = iter, cost = cost, can.include = sidx, ...)

      #--- extract sampled rows from existing ---#

      samples <- existing[outIdx, ]
    }
  } else {
    message("Sub-sampling based on 'raster' distributions.")

    #--- determine crs of input raster ---#
    crs <- terra::crs(raster, proj = TRUE)

    # #--- extract covariates data from raster ---#
    vals_r <- terra::as.data.frame(raster, xy = TRUE, row.names = FALSE) %>%
      dplyr::rename(
        X = x,
        Y = y
      ) %>%
      stats::na.omit()

    #--- select the variables existing in raster for sampling ---#
    
    vals <- vals %>%
      dplyr::select(names(vals_r))
    
    all <- dplyr::bind_rows(vals_r, vals)


    #--- sampling ---#

    sidx <- (nrow(all) - nrow(existing) + 1):nrow(all)

    if (isTRUE(details)) {
      #--- output clhs information to be supplied in 'details' list output ---#

      clhsOut <- clhs::clhs(x = all[, names(raster)], size = nSamp, iter = iter, cost = cost, can.include = sidx, simple = FALSE, ...)

      outIdx <- clhsOut$index_samples
      #--- extract sampled rows from existing ---#

      samples <- all[outIdx, ]
    } else {
      outIdx <- clhs::clhs(x = all[, names(raster)], size = nSamp, iter = iter, cost = cost, can.include = sidx, ...)

      #--- extract sampled rows from existing ---#

      samples <- all[outIdx, ]
    }

    #--- extract sampled rows from existing ---#

    samples <- samples %>%
      sf::st_as_sf(., coords = c("X", "Y"), crs = crs)
  }

  #--- write outputs if desired ---#
  write_samples(samples = samples, filename = filename, overwrite = overwrite)

  #--- details ---#

  if (isTRUE(details)) {
    out <- list(
      samples = samples,
      population = all,
      clhsOut = clhsOut,
      plotcat = if (exists("pecdfcat")) pecdfcat,
      plot = if (exists("pecdf")) pecdf
    )

    return(out)
  } else {
    return(samples)
  }
}
