#' Stratified sampling
#'
#' @description Sampling based on a stratified raster.
#'
#' @family sample functions
#'
#' @inheritParams sample_srs
#' @inheritParams calculate_allocation
#' @param sraster spatRaster. Stratification raster to be used for sampling.
#' @param nSamp Numeric. Number of desired samples. \code{existing include} and \code{force} influence this value.
#' @param existing sf or data.frame.  Existing plot network.
#' @param include Logical. If \code{TRUE} include existing plots in \code{nSamp} total.
#' @param remove Logical. If \code{TRUE} randomly remove samples from over represented strata to meet allocated sample numbers.
#' Used only when \code{existing} and \code{include} are both \code{TRUE}.
#' @param wrow Numeric. Number of row in the focal window (default is 3).
#' @param wcol Numeric. Number of columns in the focal window (default is 3).
#' @param details Logical. If \code{FALSE} (default) output is sf object of
#' stratified samples. If \code{TRUE} return a list
#' where \code{$details} additional sampling information and \code{$raster}
#' is an sf object of stratified samples.
#' @param plot Logical. Plots existing (circles) and new (crosses) samples.
#'
#' @importFrom methods is
#'
#' @return An sf object with \code{nSamp} stratified samples.
#'
#' @examples
#' #--- Load raster and access files ---#
#' r <- system.file("extdata", "kmeans.tif", package = "sgsR")
#' sr <- terra::rast(r)
#'
#' a <- system.file("extdata", "roads.shp", package = "sgsR")
#' ac <- sf::st_read(a)
#'
#' e <- system.file("extdata", "existing.shp", package = "sgsR")
#' e <- sf::st_read(e)
#'
#' #--- perform stratified sampling random sampling ---#
#' sample_strat(
#'   sraster = sr,
#'   nSamp = 200,
#'   plot = TRUE
#' )
#'
#' #--- perform stratified sampling random sampling ---#
#' sample_strat(
#'   sraster = sr,
#'   nSamp = 200,
#'   plot = TRUE,
#'   force = TRUE
#' )
#'
#' #--- extract strata values to existing samples ---#
#' e.sr <- extract_strata(sraster = sr, existing = e)
#'
#' sample_strat(
#'   sraster = sr,
#'   nSamp = 200,
#'   access = ac,
#'   existing = e.sr,
#'   mindist = 200,
#'   buff_inner = 50,
#'   buff_outer = 200
#' )
#'
#' sample_strat(
#'   sraster = sr,
#'   nSamp = 200,
#'   access = ac,
#'   buff_inner = 50,
#'   buff_outer = 200,
#'   filename = tempfile(fileext = ".shp")
#' )
#'
#' #--- Load mraster for optimal allocation ---#
#' mr <- system.file("extdata", "wall_metrics.tif", package = "sgsR")
#' mr <- terra::rast(mr)
#'
#' sample_strat(
#'   sraster = sr,
#'   nSamp = 200,
#'   allocation = "optim",
#'   mraster = mr$zq90,
#'   access = ac,
#'   buff_inner = 50,
#'   buff_outer = 200,
#'   filename = tempfile(fileext = ".shp")
#' )
#' @author Tristan R.H. Goodbody & Martin Queinnec
#'
#' @note
#' The sampling is performed in 2 stages:
#' \enumerate{
#'
#' \item \code{Rule 1} - Sample within grouped stratum pixels defined within the
#' \code{wrow, wcol} parameters
#'
#' \item \code{Rule 2} - If no samples exist to satisfy Rule 1
#'  individual stratum pixels are sampled.
#'
#'  The rule applied to allocate each sample is defined in the \code{rule} attribute of output samples.
#'
#' }
#' @references
#' Queinnec, M., White, J. C., & Coops, N. C. (2021).
#' Comparing airborne and spaceborne photon-counting LiDAR canopy
#' structural estimates across different boreal forest types.
#' Remote Sensing of Environment, 262 (August 2020), 112510.
#' https://doi.org/10.1016/j.rse.2021.112510
#'
#' @export

sample_strat <- function(sraster,
                         nSamp,
                         allocation = "prop",
                         force = FALSE,
                         mraster = NULL,
                         mindist = NULL,
                         existing = NULL,
                         include = FALSE,
                         remove = FALSE,
                         access = NULL,
                         buff_inner = NULL,
                         buff_outer = NULL,
                         wrow = 3,
                         wcol = 3,
                         plot = FALSE,
                         details = FALSE,
                         filename = NULL,
                         overwrite = FALSE) {

  #--- check for required packages ---#

  #--- both packages required only if 'mindist' is defined ---#
  if (!is.null(mindist)) {
    if (!requireNamespace("spatstat", quietly = TRUE)) {
      stop("Package \"spatstat\" needed for this function to work. Please install it.",
        call. = FALSE
      )
    }

    #--- check for required packages ---#
    if (!requireNamespace("spatstat.geom", quietly = TRUE)) {
      stop("Package \"spatstat.geom\" needed for this function to work. Please install it.",
        call. = FALSE
      )
    }
  }

  #--- Set global vars ---#
  x <- y <- cell <- NULL

  #--- Error management ---#
  if (!inherits(sraster, "SpatRaster")) {
    stop("'sraster' must be type SpatRaster", call. = FALSE)
  }

  if (any(!c("strata") %in% names(sraster))) {
    stop("'sraster must have a layer named 'strata'", call. = FALSE)
  }

  if (!is.null(mindist)) {
    if (!is.numeric(mindist)) {
      stop("'mindist' must be type numeric", call. = FALSE)
    }
  }

  if (!is.numeric(nSamp)) {
    stop("'nSamp' must be type numeric", call. = FALSE)
  }

  if (!is.logical(include)) {
    stop("'include' must be type logical", call. = FALSE)
  }

  if (!is.logical(remove)) {
    stop("'remove' must be type logical", call. = FALSE)
  }

  if (!is.logical(force)) {
    stop("'force' must be type logical", call. = FALSE)
  }

  if (!is.numeric(wrow)) {
    stop("'wrow' must be type numeric", call. = FALSE)
  }

  if (!is.numeric(wcol)) {
    stop("'wcol' must be type numeric", call. = FALSE)
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical", call. = FALSE)
  }

  if (!is.logical(details)) {
    stop("'details' must be type logical", call. = FALSE)
  }

  #--- if the sraster has multiple bands subset the band named strata ---#
  if (terra::nlyr(sraster) > 1) {
    sraster <- terra::subset(sraster, "strata")
  }

  #--- determine crs of input sraster ---#
  crs <- terra::crs(sraster, proj = TRUE)

  #--- if existing samples are provided ensure they are in the proper format ---#

  if (is.null(existing)) {
    if (isTRUE(include)) {
      stop("'existing' must be provided when 'include' == TRUE", call. = FALSE)
    }

    if (isTRUE(remove)) {
      stop("'existing' must be provided when 'remove' == TRUE", call. = FALSE)
    }

    #--- if existing samples do not exist make an empty data.frame called addSamples ---#
    addSamples <- data.frame(cell = NA, strata = NA, X = NA, Y = NA)
    extraCols <- character(0)
  } else {

    #--- existing must be either a data.frame or an sf object with columns names 'X' 'Y' 'strata' ---#

    if (!inherits(existing, "data.frame") && !inherits(existing, "sf")) {
      stop("'existing' must be a data.frame or sf object", call. = FALSE)
    }

    if (any(!c("strata") %in% names(existing))) {
      stop("'existing' must have an attribute named 'strata'. Consider using extract_strata().", call. = FALSE)
    }

    if (inherits(sf::st_geometry(existing), "sfc_POINT")) {

      #--- if existing is an sf object extract the coordinates and the strata vector ---#

      exist_xy <- sf::st_coordinates(existing)

      strata <- existing$strata

      existing <- as.data.frame(cbind(strata, exist_xy))
    } else {
      stop("'existing' geometry type must be 'sfc_POINT'", call. = FALSE)
    }

    #--- if existing samples do exist ensure proper naming convention ---#

    if (any(!c("X", "Y") %in% colnames(existing))) {

      #--- if coordinate column names are lowercase change them to uppercase to match requirements ---#

      if (any(c("x", "y") %in% colnames(existing))) {
        existing <- existing %>%
          dplyr::rename(
            X = x,
            Y = y
          )

        message("'existing' column coordinate names are lowercase - converting to uppercase")
      } else {

        #--- if no x/y columns are present stop ---#

        stop("'existing' must have columns named 'X' and 'Y'", call. = FALSE)
      }
    }

    #--- add cell value for future checking for duplicate samples ---#

    existing$cell <- NA

    addSamples <- existing
  }

  extraCols <- colnames(existing)[!colnames(existing) %in% c("cell", "X", "Y", "strata")]

  # Transform strata to numeric if factor
  if (is(addSamples$strata, "factor")) {
    addSamples$strata <- as.numeric(as.character(addSamples$strata))
  }

  #--- determine number of samples for each strata ---#

  if (isTRUE(include)) {
    message("'existing' samples being included in 'nSamp' total")

    toSample <- calculate_allocation(
      sraster = sraster,
      nSamp = nSamp,
      existing = existing,
      force = force,
      allocation = allocation,
      mraster = mraster
    )
  } else {
    toSample <- calculate_allocation(
      sraster = sraster,
      nSamp = nSamp,
      force = force,
      allocation = allocation,
      mraster = mraster
    )
  }


  #--- determine access buffers ---#

  if (!missing(access)) {

    #--- error handling in the presence of 'access' ---#
    if (!inherits(access, "sf")) {
      stop("'access' must be an 'sf' object", call. = FALSE)
    }

    if (!inherits(sf::st_geometry(access), "sfc_MULTILINESTRING") && !inherits(sf::st_geometry(access), "sfc_LINESTRING")) {
      stop("'access' geometry type must be 'LINESTRING' or 'MULTILINESTRING'", call. = FALSE)
    }

    if (buff_inner > buff_outer) {
      stop("'buff_inner' must be < 'buff_outer'", call. = FALSE)
    }

    access_buff <- mask_access(raster = sraster, access = access, buff_inner = buff_inner, buff_outer = buff_outer)

    raster_masked <- access_buff$rast
  }

  ####################################
  #--- Start of sampling function ---#
  ####################################

  for (i in 1:nrow(toSample)) {
    s <- as.numeric(toSample[i, 1])
    n <- as.numeric(toSample[i, 2])

    message(glue::glue("Processing strata : {s}"))

    #--- if the number of samples required is equal to zero (if `include == TRUE`) just keep existing samples only ---#
    if (n == 0) {

      #--- Initiate number of sampled cells ---#
      add_strata <- addSamples %>%
        dplyr::filter(strata == s)

      if (nrow(add_strata) > 0) {
        add_strata$type <- "existing"

        if (!"rule" %in% colnames(add_strata)) {
          add_strata$rule <- "existing"
        }
      }

      message(glue::glue("Strata : {s} required no sample additions. Keeping all existing samples."))
    } else if (n > 0) {
      #--- mask for individual strata ---#

      strata_m <- terra::mask(sraster,
        mask = sraster,
        maskvalues = s,
        inverse = TRUE
      )
      names(strata_m) <- "strata"

      #--- if access line polygon is specified create inner and outer buffers

      if (!missing(access)) {
        strata_m_buff <- terra::mask(strata_m,
          mask = access_buff$buff
        )

        sampAvail <- sum(!is.na(terra::values(strata_m_buff)))

        if (sampAvail > n) {
          message(
            glue::glue("Buffered area contains {sampAvail} available candidates. Sampling to reach {n} samples starting.")
          )

          #--- rename to original strata sraster that will be used for sampling ---#
          strata_m <- strata_m_buff

          #--- if there are no samples to take within the specified 'buff_outer' distance extend buffer until values are found ---#
        } else {
          stop("Insufficient candidate samples within the buffered access extent. Consider altering buffer widths.", call. = FALSE)
        }
      }

      ### --- sampling ---###

      ### --- RULE 1: select only cells surrounded by cells with same strata ---###

      #--- Define focal window ---#
      w <- matrix(1 / (wrow * wcol), wrow, wcol)

      suppressWarnings(strata_m_clust <-
        terra::focal(
          strata_m,
          w = w,
          na.rm = FALSE
        ))
      names(strata_m_clust) <- "strata"

      #--- Initiate number of sampled cells ---#
      add_strata <- addSamples %>%
        dplyr::filter(strata == s)

      if (nrow(add_strata) > 0) {
        add_strata$type <- "existing"

        if (!"rule" %in% colnames(add_strata)) {
          add_strata$rule <- "existing"
        }
      }

      #--- create indices for all, NA, and valid sampling candidates ---#

      idx_all <- 1:terra::ncell(strata_m_clust)
      idx_na <- !complete.cases(terra::values(strata_m_clust))
      validCandidates <- idx_all[!idx_na]

      #--- Rule 1 sampling ---#
      nCount <- 0 # Number of sampled cells

      # While loop for RULE 1
      while (length(validCandidates) > 0 & nCount < n) {
        #-- identify potential sample from candidates ---#
        smp <- sample(1:length(validCandidates), size = 1)

        smp_cell <- validCandidates[smp]

        #--- Remove sampled cell from validCandidates so that it cannot be sampled again later ---#
        validCandidates <- validCandidates[-smp]

        #--- extract coordinates and sample details ---#

        add_temp <- data.frame(
          cell = smp_cell,
          X = terra::xFromCell(strata_m_clust, smp_cell),
          Y = terra::yFromCell(strata_m_clust, smp_cell),
          strata = strata_m_clust[smp_cell]
        )

        #--- populate add_temp with values ---#
        add_temp$type <- "new"
        add_temp$rule <- "rule1"
        add_temp[, extraCols] <- NA

        #--- If add_strata is empty, sampled cell accepted ---#

        if (nrow(add_strata) == 0) {
          add_strata <- add_temp[, c("cell", "X", "Y", "strata", "type", "rule", extraCols)]

          nCount <- nCount + 1

          #--- If add_strata isn't empty, check distance with all other sampled cells in strata ---#
        }

        if (!is.null(mindist)) {
          dist <- spatstat.geom::crossdist(add_temp$X, add_temp$Y, add_strata$X, add_strata$Y)

          #--- If all less than 'mindist' - accept sampled cell otherwise reject ---#
          if (all(as.numeric(dist) > mindist)) {
            add_strata <- rbind(add_strata, add_temp[, c("cell", "X", "Y", "strata", "type", "rule", extraCols)])

            nCount <- nCount + 1
          }
        } else {

          #--- if mindist is not defined ---#

          if (add_temp$cell %in% add_strata$cell) next

          add_strata <- rbind(add_strata, add_temp[, c("cell", "X", "Y", "strata", "type", "rule", extraCols)])

          nCount <- nCount + 1
        }
      }

      ### --- RULE 2 sampling ---###

      if (nCount < n) {
        idx_all <- 1:terra::ncell(strata_m)
        idx_na <- !complete.cases(terra::values(strata_m))
        validCandidates <- idx_all[!idx_na]

        while (length(validCandidates) > 0 & nCount < n) {

          #-- identify potential sample from candidates ---#
          smp <- sample(1:length(validCandidates), size = 1)

          smp_cell <- validCandidates[smp]

          #--- Remove sampled cell from validCandidates so that it cannot be sampled again later ---#

          validCandidates <- validCandidates[-smp]

          #--- extract coordinates and sample details ---#

          add_temp <- data.frame(
            cell = smp_cell,
            X = terra::xFromCell(strata_m, smp_cell),
            Y = terra::yFromCell(strata_m, smp_cell),
            strata = validCandidates[smp_cell]
          )

          add_temp$rule <- "rule2"
          add_temp$type <- "new"
          add_temp[, extraCols] <- NA
          add_temp$strata <- s

          if (nrow(add_strata) == 0) {
            add_strata <- add_temp[, c("cell", "X", "Y", "strata", "type", "rule", extraCols)]

            nCount <- nCount + 1

            #--- If add_strata isn't empty, check distance with all other sampled cells in strata ---#
          }

          if (!is.null(mindist)) {
            dist <- spatstat.geom::crossdist(add_temp$X, add_temp$Y, add_strata$X, add_strata$Y)

            #--- If all less than 'mindist' - accept sampled cell otherwise reject ---#
            if (all(as.numeric(dist) > mindist)) {
              add_strata <- rbind(add_strata, add_temp[, c("cell", "X", "Y", "strata", "type", "rule", extraCols)])

              nCount <- nCount + 1
            }
          } else {

            #--- if mindist is not defined ---#

            if (add_temp$cell %in% add_strata$cell) next

            add_strata <- rbind(add_strata, add_temp[, c("cell", "X", "Y", "strata", "type", "rule", extraCols)])

            nCount <- nCount + 1
          }
        }
      }

      if (nCount < n) {
        message(sprintf("Strata %s: couldn't select required number of samples: %i instead of %i \n", s, nCount, n))
      }

      #--- if number of samples is < 0 based on `include` parameter ---#
    } else if (n < 0) {
      if (isTRUE(remove)) {

        #--- need to remove samples from over represented strata ---#

        #--- sample total needed from existing ---#
        need <- as.numeric(toSample[i, 3])

        message(glue::glue("'include = TRUE & remove = TRUE' - Stratum {s} overrepresented - {abs(n)} samples removed."))

        add_strata <- addSamples %>%
          dplyr::filter(strata == s) %>%
          dplyr::sample_n(need)

        #--- add type and rule attributes ---#

        add_strata$type <- "existing"
        add_strata$rule <- "existing"
      } else {
        message(glue::glue("'include = TRUE & remove = FALSE' - Stratum {s} overrepresented by {abs(n)} samples but have not been removed. Expect a higher total 'nSamp' in output"))
        #--- keep over represented samples in dataset ---#
        add_strata <- addSamples %>%
          dplyr::filter(strata == s)

        if (nrow(add_strata) > 0) {
          add_strata$type <- "existing"

          if (!"rule" %in% colnames(add_strata)) {
            add_strata$rule <- "existing"
          }
        }
      }
    }

    # Create out object if first iteration of loop
    # Else just rbind output with what has been processed in the loop
    if (i == 1) {
      out <- add_strata
    } else {
      out <- rbind(out, add_strata)
    }
  }

  #--- convert coordinates to a spatial points object ---#
  samples <- out %>%
    dplyr::select(-cell) %>%
    as.data.frame() %>%
    sf::st_as_sf(., coords = c("X", "Y"))

  #--- assign sraster crs to spatial points object ---#
  sf::st_crs(samples) <- crs

  #--- plot the raster and samples if desired ---#

  if (isTRUE(plot)) {

    #--- if existing is not provided plot the masked raster ---#

    if (missing(existing)) {

      #--- if access is also missing plot the full sraster extent ---#

      if (missing(access)) {
        terra::plot(sraster[[1]])
        suppressWarnings(terra::plot(samples, add = T, col = "black"))

        #--- if access is provided plot the masked access sraster ---#
      } else {
        terra::plot(sraster[[1]])
        suppressWarnings(terra::plot(access_buff$buff, add = T, border = c("gray30"), col = "gray10", alpha = 0.1))
        suppressWarnings(terra::plot(samples, add = T, col = "black"))
      }

      #--- if existing is provided plot the full raster ---#
    } else {

      #--- plot input sraster and random samples ---#

      terra::plot(sraster[[1]])
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

    output <- list(sampleDist = toSample, samples = samples)

    #--- output samples dataframe ---#
    return(output)
  } else {

    #--- just output raster ---#

    return(samples)
  }
}
