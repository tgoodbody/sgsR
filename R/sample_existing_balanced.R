#' Sample Existing Data Using Balanced Sampling
#'
#' This function samples a given set of existing data using balanced sampling techniques,
#' which ensures that each stratum or subgroup of data is proportionally represented in the sample.
#'
#' @inheritParams sample_balanced
#' @inheritParams extract_strata 
#' @param ... Additional arguments to pass to the selected sampling algorithm. 
#' This is leveraged when used by \code{sample_existing()} internally
#' @return An \code{sf} object that is a sub-sample of \code{existing} 
#' @keywords internal
sample_existing_balanced <- function(existing,
                                     nSamp,
                                     algorithm = "lpm2_kdtree",
                                     p = NULL,
                                     filename = NULL,
                                     overwrite = FALSE,
                                     ...) {
  #--- Set global vars ---#
  x <- y <- X <- Y <- strata <- NULL

  if (!is.character(algorithm)) {
    stop("'algorithm' must be type character.", call. = FALSE)
  }

  #--- list all available algorithms to determine if a valid one has been supplied ---#
  algs <- c("lpm2_kdtree", "lcube", "lcubestratified")

  if (!algorithm %in% algs) {
    stop("Unknown algorithm specified. Please use one of 'lpm2_kdtree', 'lcube', 'lcubestratified'.", call. = FALSE)
  }

  #--- get existing values ---#
  vals <- coords_existing(existing)

  vals_m <- as.matrix(vals %>% dplyr::select(-X, -Y))

  N <- nrow(vals)

  if (is.null(p)) {
    p <- rep(nSamp / N, N)
  } else {
    if (!is.numeric(p)) {
      stop("'p' must be type numeric.", call. = FALSE)
    }
    if (length(p) != N) {
      stop(paste0("'p' have a length of ", N, "."), call. = FALSE)
    }
  }
  if (algorithm == "lpm2_kdtree") {
    sampled <- SamplingBigData::lpm2_kdtree(prob = p, x = vals_m)
  }
  if (algorithm == "lcube") {
    sampled <- BalancedSampling::lcube(
      prob = p, Xspread = vals_m,
      Xbal = cbind(p)
    )
  }
  if (algorithm == "lcubestratified") {
    if (!"strata" %in% names(existing)) {
      stop("'existing' must have a variable named 'strata' to use the 'lcubestratified' algorithm.",
        call. = FALSE
      )
    }
    strata_v <- as.vector(vals$strata)
    vals_m <- as.matrix(dplyr::select(vals, -X, -Y, -strata))
    sampled <- BalancedSampling::lcubestratified(
      prob = p,
      Xspread = vals_m, Xbal = cbind(p), integerStrata = strata_v
    )
    sampled <- (1:N)[sampled == 1]
  }
  samples <- vals[sampled, ]
  samples <- samples %>%
    as.data.frame() %>%
    sf::st_as_sf(.,
      coords = c("X", "Y"),
      crs = sf::st_crs(existing)
    )

  if (!is.null(filename)) {
    if (!is.character(filename)) {
      stop("'filename' must be a file path character string.", call. = FALSE)
    }

    if (!is.logical(overwrite)) {
      stop("'overwrite' must be type logical.", call. = FALSE)
    }

    if (file.exists(filename) & isFALSE(overwrite)) {
      stop(paste0("'", filename, "' already exists and overwrite = FALSE."), call. = FALSE)
    }

    sf::st_write(samples, filename, delete_layer = overwrite)
    message("Output samples written to disc.")
  }

  return(samples)
}
