#' Stratified sampling
#'
#' @description Sampling based on a stratified raster.
#'
#' @family sample functions
#'
#' @inheritParams sample_srs
#' @inheritParams calculate_allocation
#' @param method Character. Sampling design approach \code{"Queinnec"} (default) or \code{"random"}. \code{"Queinnec"} method is
#' described in notes below. \code{"random"} performs traditional stratified random sampling where probability to sample each
#' cell within each stratum is equal assuming default parameters for \code{mindist}. \code{existing, include, remove} are ignored when \code{method = "random"}.
#' @param sraster spatRaster. Stratification raster to be used for sampling.
#' @param nSamp Numeric. Number of desired samples. \code{existing}, \code{include} and \code{force} influence this value.
#' @param existing sf 'POINT' or data.frame.  Existing plot network.
#' @param include Logical. If \code{TRUE} include \code{existing} plots in \code{nSamp} total.
#' @param remove Logical. If \code{TRUE} randomly remove samples from over represented strata to meet allocated sample numbers.
#' Used only when \code{existing} and \code{include} are both \code{TRUE}.
#' @param wrow Numeric. Number of row in the focal window (\code{default = 3}).
#' @param wcol Numeric. Number of columns in the focal window (\code{default = 3}).
#' @param details Logical. If \code{FALSE} (default) output is sf object of
#' stratified samples. If \code{TRUE} return a list
#' where \code{$details} additional sampling information and \code{$raster}
#' is an sf object of stratified samples.
#' @param plot Logical. Plots samples of type `existing` (if provided; croses) and `new` (circles) along with \code{sraster}.
#'
#' @importFrom methods is
#'
#' @return An sf object with \code{nSamp} stratified samples.
#'
#' @examples
#' #--- Load raster and access files ---#
#' r <- system.file("extdata", "sraster.tif", package = "sgsR")
#' sr <- terra::rast(r)
#'
#'
#' #--- perform stratified sampling random sampling ---#
#' sraster <- sample_strat(
#'   sraster = sr,
#'   nSamp = 50
#' )
#'
#' @author Tristan R.H. Goodbody & Martin Queinnec
#'
#' @note
#' The sampling is performed in 2 stages when \code{method = "Queinnec"}:
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
#'
#' \code{existing} may contain samples that fall in \code{sraster} cells that are `NA`. If this occurs and \code{include = TRUE}, `NA` samples
#' are separated during sampling and re-appended at the end of the sampling process.
#'
#' If the \code{sraster} provided contains factor values, the algorithm will automatically convert these into the numeric factor levels and
#' perform sampling using those values. The categories (factor values) will be extracted and appended to the algorithm output as the `category` attribute.
#'
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
                         method = "Queinnec",
                         weights = NULL,
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
  #--- Set global vars ---#
  x <- y <- cell <- cats <- rule <- NULL

  #--- Error management ---#
  if (!inherits(sraster, "SpatRaster")) {
    stop("'sraster' must be type SpatRaster.", call. = FALSE)
  }

  if (!is.numeric(nSamp)) {
    stop("'nSamp' must be type numeric.", call. = FALSE)
  }

  if (!is.character(method)) {
    stop("'method' must be type character.", call. = FALSE)
  }

  if (!method %in% c("random", "Queinnec")) {
    stop("'method' must be one of 'random' or 'Queinnec'.", call. = FALSE)
  }

  if (!is.null(mindist)) {
    if (!is.numeric(mindist)) {
      stop("'mindist' must be type numeric.", call. = FALSE)
    }
  }

  if (!is.logical(include)) {
    stop("'include' must be type logical.", call. = FALSE)
  }

  if (!is.logical(remove)) {
    stop("'remove' must be type logical.", call. = FALSE)
  }

  if (!is.logical(force)) {
    stop("'force' must be type logical.", call. = FALSE)
  }

  if (!is.numeric(wrow)) {
    stop("'wrow' must be type numeric.", call. = FALSE)
  }

  if (!is.numeric(wcol)) {
    stop("'wcol' must be type numeric.", call. = FALSE)
  }

  if ((wrow %% 2) == 0) {
    stop("'wrow' must be an odd number.", call. = FALSE)
  }

  if ((wcol %% 2) == 0) {
    stop("'wcol' must be an odd number.", call. = FALSE)
  }

  if (!is.logical(plot)) {
    stop("'plot' must be type logical.", call. = FALSE)
  }

  if (!is.logical(details)) {
    stop("'details' must be type logical.", call. = FALSE)
  }

  #--- check if `sraster` contains factor values and if so generate its category list to amend later ---#

  if (!is.null(terra::cats(sraster)[[1]])) {
    message("'sraster' has factor values. Converting to allow mapping.")

    #--- change suggested by R Hijmans ---#
    sraster_cats <- cats(sraster) %>%
      as.data.frame()
    colnames(sraster_cats)[1] <- "value"
  }

  #--- determine crs of input sraster ---#
  crs <- terra::crs(sraster)

  if (method == "Queinnec") {
    message("Using 'Queinnec' sampling method.")

    #--- if existing samples are provided ensure they are in the proper format ---#

    if (is.null(existing)) {
      if (isTRUE(include)) {
        stop("'existing' must be provided when 'include = TRUE'.", call. = FALSE)
      }

      if (isTRUE(remove)) {
        stop("'existing' must be provided when 'remove = TRUE'.", call. = FALSE)
      }

      #--- if existing samples do not exist make an empty data.frame called addSamples ---#
      addSamples <- data.frame(cell = NA, strata = NA, X = NA, Y = NA)
      extraCols <- character(0)
    } else {
      #--- existing must be either a data.frame or an sf object with columns names 'X' 'Y' 'strata' ---#

      if (!inherits(existing, "data.frame") && !inherits(existing, "sf")) {
        stop("'existing' must be a data.frame or sf object.", call. = FALSE)
      }

      if (inherits(existing, "sf")) {
        if (inherits(sf::st_geometry(existing), "sfc_POINT")) {
          #--- determine crs of existing ---#
          crs <- sf::st_crs(existing)

          #--- if existing is an sf object extract the coordinates and the strata vector ---#

          exist_xy <- sf::st_coordinates(existing)

          strata <- existing$strata

          existing <- as.data.frame(cbind(strata, exist_xy))
        } else {
          stop("'existing' geometry type must be 'sfc_POINT'.", call. = FALSE)
        }
      }

      if (any(!c("strata") %in% names(existing))) {
        stop("'existing' must have an attribute named 'strata'. Consider using extract_strata().", call. = FALSE)
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

          message("'existing' column coordinate names are lowercase - converting to uppercase.")
        } else {
          #--- if no x/y columns are present stop ---#

          stop("'existing' must have columns named 'X' and 'Y'.", call. = FALSE)
        }
      }

      #--- add cell value for future checking for duplicate samples ---#

      existing$cell <- NA

      addSamples <- existing
    }


    extraCols <- colnames(existing)[!colnames(existing) %in% c("cell", "X", "Y", "strata")]
  } else {
    message("Using 'random' sampling method. Ignoring 'existing', 'include', 'remove' if provided.")
    
    addSamples <- data.frame(cell = NA, strata = NA, X = NA, Y = NA)
    extraCols <- character(0)

    existing <- NULL
    include <- NULL
    remove <- NULL
  }

  #--- determine number of samples for each strata ---#

  if (isTRUE(include)) {
    message("'existing' samples being included in 'nSamp' total.")

    toSample <- calculate_allocation(
      sraster = sraster,
      nSamp = nSamp,
      weights = weights,
      existing = existing,
      force = force,
      allocation = allocation,
      mraster = mraster
    )
  } else {
    toSample <- calculate_allocation(
      sraster = sraster,
      nSamp = nSamp,
      weights = weights,
      force = force,
      allocation = allocation,
      mraster = mraster
    )
  }

  #--- determine access buffers ---#

  if (!missing(access)) {
    access_buff <- mask_access(raster = sraster, access = access, buff_inner = buff_inner, buff_outer = buff_outer)

    raster_masked <- access_buff$rast
  }

  #--- Define focal window ---#
  w <- matrix(1 / (wrow * wcol), wrow, wcol)

  ####################################
  #--- Start of sampling function ---#
  ####################################

  for (i in 1:nrow(toSample)) {
    s <- as.numeric(toSample[i, 1])
    n <- as.numeric(toSample[i, 2])

    message(paste0("Processing strata : ", s))

    #--- use stratified RANDOM sampling or "Queinnec" method ---#

    if (method == "random") {
      if (n == 0) {
        message("No samples needed.")
      } else if (n > 0) {
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
              paste0("Buffered area contains ", sampAvail, " available candidates. Sampling to reach ", n, " starting.")
            )

            #--- rename to original strata sraster that will be used for sampling ---#
            strata_m <- strata_m_buff

            #--- if there are no samples to take within the specified 'buff_outer' distance extend buffer until values are found ---#
          } else {
            stop("Insufficient candidate samples within the buffered access extent. Consider altering buffer widths.", call. = FALSE)
          }
        }
        
        #--- initiate number of sampled cells ---#
        add_strata <- addSamples %>%
          dplyr::filter(strata == s)

        #--- ensure that sample units from previous strata are appended for distance checking ---#
        if(!is.null(mindist)){
          if(exists("out")){
            
            add_strata <- rbind(out, add_strata)
            
          }
        }
        
        add_strata <- strat_rule2(
          n = n,
          s = s,
          add_strata = add_strata,
          nCount = 0,
          strata_m = strata_m,
          extraCols = extraCols,
          mindist = mindist
        ) %>%
          filter(strata == s)
        
      }
    }

    if (method == "Queinnec") {
      #--- if the number of samples required is equal to zero (if `include = TRUE`) just keep existing samples only ---#
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

        message(paste0("Strata : ", s, " required no sample additions. Keeping all existing samples."))
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
              paste0("Buffered area contains ", sampAvail, " available candidates. Sampling to reach ", n, " starting.")
            )

            #--- rename to original strata sraster that will be used for sampling ---#
            strata_m <- strata_m_buff

            #--- if there are no samples to take within the specified 'buff_outer' distance extend buffer until values are found ---#
          } else {
            stop("Insufficient candidate samples within the buffered access extent. Consider altering buffer widths.", call. = FALSE)
          }
        }

        ### --- sampling ---###

        suppressWarnings(strat_mask <-
          terra::focal(
            strata_m,
            w = w,
            na.rm = FALSE
          ))
        names(strat_mask) <- "strata"

        #--- Initiate number of sampled cells ---#
        add_strata <- addSamples %>%
          dplyr::filter(strata == s)

        if (nrow(add_strata) > 0) {
          add_strata$type <- "existing"

          if (!"rule" %in% colnames(add_strata)) {
            add_strata$rule <- "existing"
          }
        }
        
        #--- ensure that sample units from previous strata are appended for distance checking ---#
        if(!is.null(mindist)){
          if(exists("out")){
            
            add_strata <- rbind(out, add_strata)
            
          }
        }

        #--- Rule 1 sampling ---#

        r1 <- strat_rule1(
          n = n,
          i = i,
          s = s,
          strat_mask = strat_mask,
          add_strata = add_strata,
          extraCols = extraCols,
          mindist = mindist
        )

        #--- Rule 2 sampling ---#

        add_strata <- strat_rule2(
          n = n,
          s = s,
          add_strata = r1$add_strata,
          nCount = r1$nCount,
          strata_m = strata_m,
          extraCols = extraCols,
          mindist = mindist
        ) %>%
          filter(strata == s)
          

        #--- if number of samples is < 0 based on `include` parameter ---#
      } else if (n < 0) {
        if (isTRUE(remove)) {
          #--- need to remove samples from over represented strata ---#

          #--- sample total needed from existing ---#
          need <- as.numeric(toSample[i, 3])

          message(paste0("'include = TRUE & remove = TRUE' - Stratum ", s, " overrepresented - ", abs(n), " samples removed."))

          add_strata <- addSamples %>%
            dplyr::filter(strata == s) %>%
            dplyr::sample_n(need)

          #--- add type and rule attributes ---#

          add_strata$type <- "existing"
          add_strata$rule <- "existing"
        } else {
          message(paste0("'include = TRUE & remove = FALSE' - Stratum ", s, " overrepresented by ", abs(n), " samples but have not been removed. Expect a higher total 'nSamp' in output."))
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
    }

    # Create out object if first iteration of loop
    # Else just rbind output with what has been processed in the loop
    if (i == 1) {
      out <- add_strata
    } else {
      out <- rbind(out, add_strata)
    }
  }

  #--- check if samples fall in areas where stratum values are NA ---#
  if (!is.null(existing)) {
    if (any(!complete.cases(existing$strata))) {
      na_only <- existing %>%
        dplyr::filter(!complete.cases(strata)) %>%
        dplyr::select(-cell)

      samples_NA <- na_only %>%
        dplyr::mutate(
          type = "existing",
          rule = NA
        )

      #--- convert coordinates to a spatial points object ---#
      samples <- out %>%
        dplyr::select(-cell) %>%
        as.data.frame() %>%
        rbind(., samples_NA) %>%
        sf::st_as_sf(., coords = c("X", "Y"), crs = crs)
    } else {
      #--- convert coordinates to a spatial points object ---#
      samples <- out %>%
        dplyr::select(-cell) %>%
        as.data.frame() %>%
        sf::st_as_sf(., coords = c("X", "Y"), crs = crs)
    }
  } else {
    #--- convert coordinates to a spatial points object ---#
    samples <- out %>%
      dplyr::select(-cell) %>%
      as.data.frame() %>%
      sf::st_as_sf(., coords = c("X", "Y"), crs = crs)
    
    if(method == "random"){
      
      samples <- samples %>%
        dplyr::select(-rule)
      
    }
  
  }

  #--- plot the raster and samples if desired ---#

  if (isTRUE(plot)) {
    if (method == "random") {
      if (missing(access)) {
        terra::plot(sraster[[1]])
        suppressWarnings(terra::plot(samples, add = T, col = "black"))

        #--- if access is provided plot the masked access sraster ---#
      } else {
        terra::plot(sraster[[1]])
        suppressWarnings(terra::plot(access_buff$buff, add = T, border = c("gray30"), col = "gray10", alpha = 0.1))
        suppressWarnings(terra::plot(samples, add = T, col = "black"))
      }
    } else {
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
        suppressWarnings(terra::plot(samples, add = T, col = "black", pch = ifelse(samples$type == "existing", 3, 1)))
      }
    }
  }

  if (exists("sraster_cats")) {
    #--- match label to value from categorical raster ---#
    
    samples$category <- sraster_cats$label[match(samples$strata,sraster_cats$value)]

  }

  #--- write outputs if desired ---#
  write_samples(samples = samples, filename = filename, overwrite = overwrite)

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
