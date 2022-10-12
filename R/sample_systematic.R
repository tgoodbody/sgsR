#' Systematic sampling
#'
#' @description Systematic sampling with random start point and rotation within a square or hexagonal tessellation.
#'
#' @family sample functions
#'
#' @param raster spatRaster. Raster used to define extent of fishnet grid.
#' @param cellsize Numeric. Desired cellsize for tessellation.
#' @param square Logical. Tessellation shape. Default is regular square grid,
#' if \code{FALSE} hexagons are used.
#' @param location Character. Sample location within tessellation. \code{Default = "centers"})
#' returns samples at tessellation centers, \code{"corners"} - corners of tessellation are returned,
#' \code{"random"} - samples are randomly located within tessellations.
#' @param force Logical. Only applies when \code{location = "random"}. If \code{TRUE}, random samples are
#' forced to fall in areas where \code{raster} does not have \code{NA} values. This will considerably slow processing.
#' @param access sf. Road access network - must be lines.
#' @param buff_inner Numeric. Inner buffer boundary specifying distance
#'  from access where plots cannot be sampled.
#' @param buff_outer Numeric. Outer buffer boundary specifying distance
#'  from access where plots can be sampled.
#' @param plot Logical. Plots output strata raster with samples.
#' @param filename Character. Path to write output samples.
#' @param overwrite Logical. Choice to overwrite existing \code{filename} if it exists.
#' @param details Logical. If \code{FALSE} (default) output is sf object of
#' systematic samples. If \code{TRUE} returns a list of sf objects where \code{tessellation}
#' is the tessellation grid for sampling, and \code{samples} are the systematic samples.
#' @param ... Additional arguments for \code{\link[sf]{st_make_grid}}. Options include \code{offset} 
#' to offset grid by providing lower left coordinates.
#'
#' @return An sf object with sampled points over a tessellation.
#'
#' @note Specifying \code{location = "random"} can result in tessellations with no samples.
#' This results from \code{raster} have \code{NA} values at the random location chosen.
#' Using \code{force = TRUE} removes areas of \code{NA} from sampling entirely, but
#' considerably slows processing speeds. 
#'
#' @examples
#' #--- Load raster and access files ---#
#' r <- system.file("extdata", "sraster.tif", package = "sgsR")
#' sr <- terra::rast(r)
#'
#' a <- system.file("extdata", "access.shp", package = "sgsR")
#' ac <- sf::st_read(a)
#'
#' e <- system.file("extdata", "existing.shp", package = "sgsR")
#' e <- sf::st_read(e)
#'
#' #--- perform grid sampling ---#
#' sample_systematic(
#'   raster = sr,
#'   cellsize = 1000
#' )
#'
#' sample_systematic(
#'   raster = sr,
#'   cellsize = 1000,
#'   square = FALSE,
#'   location = "corners",
#'   plot = TRUE
#' )
#'
#' sample_systematic(
#'   raster = sr,
#'   cellsize = 1000,
#'   square = FALSE,
#'   location = "random"
#' )
#' @author Tristan R.H. Goodbody, Lukas Winiwarter
#'
#' @export


sample_systematic <- function(raster,
                              cellsize,
                              square = TRUE,
                              location = "centers",
                              force = FALSE,
                              access = NULL,
                              buff_inner = NULL,
                              buff_outer = NULL,
                              plot = FALSE,
                              filename = NULL,
                              overwrite = FALSE,
                              details = FALSE,
                              ...) {
  
  #--- Set global vars ---#
  
  ext <- geometry <- x <- pass <- NULL
  
  if (!inherits(raster, "SpatRaster")) {
    stop("'raster' must be type SpatRaster.", call. = FALSE)
  }
  
  if (!is.numeric(cellsize)) {
    stop("'cellsize' must be type numeric.", call. = FALSE)
  }
  
  if (cellsize < 0) {
    stop("'cellsize' must be > 0.", call. = FALSE)
  }
  
  if (!is.logical(plot)) {
    stop("'plot' must be type logical.", call. = FALSE)
  }
  
  if (!is.logical(square)) {
    stop("'square' must be type logical.", call. = FALSE)
  }
  
  if (!is.character(location)) {
    stop("'location' must be type character.", call. = FALSE)
  } else {
    if (!any(c("centers", "corners", "random") %in% location)) {
      stop("'location' must be one of 'centers', 'corners', or 'random'.", call. = FALSE)
    }
  }
  
  #--- determine crs of input raster ---#
  crs <- terra::crs(raster, proj = TRUE)
  
  #--- set mraster for plotting who area in case of masking ---#
  
  rasterP <- raster
  
  if (!is.null(access)) {
    
    access_buff <- mask_access(raster = raster, access = access, buff_inner = buff_inner, buff_outer = buff_outer)
    
    raster <- access_buff$rast
  }
  
  #--- convert raster extent into a polygon ---#
  
  rasterext <- sf::st_as_sf(terra::as.polygons(terra::ext(raster), crs = terra::crs(raster)))

  r <- raster
  
  res <- terra::res(x = r)[1]
  
  #--- add randomness to grid lower left coordinate locations ---#
  
  xminc <- as.numeric(terra::ext(r)[1]) + (res * sample(-100:0, 1))
  yminc <- as.numeric(terra::ext(r)[3]) + (res * sample(-100:0, 1))
  
  xx <- c(xminc,as.numeric(terra::ext(r)[2]))
  yy <- c(yminc,as.numeric(terra::ext(r)[4]))
  
  cc <- data.frame(X = xx, Y = yy)
  
  pol = sf::st_polygon(
    list(
      cbind(
        cc$X[c(1,2,2,1,1)], 
        cc$Y[c(1,1,2,2,1)])
    )
  )

  sfObj <- sf::st_sfc(pol, crs = terra::crs(r))
  
  #--- extract centroid of object for random rotation later ---#
  
  cent <- sf::st_centroid(sfObj)
  
  #--- random rotation value ---#
  
  randRot <- runif(1,0,360)
  
  #--- create tessellation ---#
  
  grid <- sf::st_make_grid(x = tran(sfObj, -randRot, cent), 
                           cellsize = cellsize, 
                           square = square, 
                           what = "polygons", 
                           crs = terra::crs(raster), 
                           ...)
  
  #--- rotate grid for plotting---#
  
  gridR <- sf::st_as_sf(tran(grid, randRot, cent),crs = terra::crs(raster))
  
  if(isTRUE(force)){
    
    if(location != "random"){
      stop("'location' must be 'random' when 'force = TRUE'", call. = FALSE)
    }
    
    message("Forcing samples to fall in non NA locations.")
    
    #--- force "random" samples to not fall in areas of no data ---#
    
    #--- generate NA raster mask ---#
    m <- terra::setValues(raster[[1]], NA)
    m[is.na(raster[[1]])] <- 1
    
    #--- convert non-NA areas to polygon ---#
    p <- sf::st_as_sf(terra::as.polygons(m)) %>%
      dplyr::select(geometry) %>%
      sf::st_combine()

    #--- if there is overlapping or slivers fix the issues ---#
    if(isFALSE(sf::st_is_valid(p))){
      
      p <- p %>%
        sf::st_make_valid()
      
    }

    p_diff <- suppressWarnings(sf::st_difference(gridR,p)) %>%
      sf::st_intersection(.,sfObj) %>%
      dplyr::filter(sf::st_is(., c("POLYGON","MULTIPOLYGON","GEOMETRYCOLLECTION")))
    
    #--- sampling --#
    samples <- sf::st_sample(p_diff, size=c(1,1), type = "random") %>%
      sf::st_sf(.)
    
    #--- check to make sure that samples intersect raster extent (cellsize check) ---#
    if(isFALSE(lengths(suppressMessages(sf::st_intersects(samples, rasterext))) > 0)){
      stop("No samples intersect with 'raster'. Ensure 'cellsize' makes sense.", call. = FALSE)
    }
    
    samples <- samples %>% 
      extract_metrics(mraster = raster[[1]], existing = ., quiet = TRUE) %>%
      stats::na.omit() %>%
      dplyr::select(geometry)
    
  } else {
    
  #--- random sampling within tesselations ---#
  
    if (location == "random") {
      location <- "centers"
      
      #--- maximum number of samples that can be selected ---#
      
      gridn <- nrow(gridR)
      
      #--- determine maximum distance sample can be moved in X / Y to remain within each tessellation ---#
      
      radius <- cellsize / 2
      
      if (square == TRUE) {

        vals <- data.frame(X = runif(gridn, -radius, radius),
                           Y = runif(gridn, -radius, radius))
      } else {
        tests <- gridn * 100
        
        X <- runif(tests, -1, 1)
        Y <- runif(tests, -.866, .866)
        
        xy <- data.frame(X = X, Y = Y)
        
        #--- test whether the sample will be within the bounds of the hexagon ---#
        xy$pass <- abs(xy$X) < 1 - .5 * (abs(xy$Y) / .866)
        
        #--- filter only values with TRUE in $pass ---#
        vals <- xy %>%
          dplyr::filter(pass == TRUE) %>%
          dplyr::slice_sample(., n = nrow(gridR)) %>%
          dplyr::mutate(
            X = X * radius,
            Y = Y * radius
          )
      }
    
    #--- create grid and locate samples ---#
    
    samples <- sf::st_centroid(gridR) %>%
      dplyr::rename(geometry = x) %>%
      sf::st_coordinates(.) %>%
      as.data.frame() %>%
      #--- apply random movement by row ---#
      dplyr::mutate(
        X = (vals$X * cos((randRot / 180) * pi)) + (vals$Y * sin((randRot / 180) * pi)) + X,
        Y = (vals$X * -sin((randRot / 180) * pi)) + (vals$Y * cos((randRot / 180) * pi)) + Y 
      ) %>%
      sf::st_as_sf(., coords = c("X", "Y"),
                   crs = terra::crs(raster))
      
    #--- check to make sure that samples intersect raster extent (cellsize check) ---#
    if(isFALSE(lengths(suppressMessages(sf::st_intersects(samples, rasterext))) > 0)){
      stop("No samples intersect with 'raster'. Ensure 'cellsize' makes sense.", call. = FALSE)
    }
      
      #--- need to extract a metric to determine if values are NA ---#
      samples <- samples %>% 
        extract_metrics(mraster = raster[[1]], existing = ., quiet = TRUE) %>%
      #--- remove samples with NA ---#
        stats::na.omit() %>%
        dplyr::select(geometry)
    
  } else if (location == "corners"){
    
    samples <- sf::st_cast(gridR, "POINT")%>%
      dplyr::rename(geometry = x)
    
    #--- check to make sure that samples intersect raster extent (cellsize check) ---#
    if(isFALSE(lengths(suppressMessages(sf::st_intersects(samples, rasterext))) > 0)){
      stop("No samples intersect with 'raster'. Ensure 'cellsize' makes sense.", call. = FALSE)
    }
      
    samples <- samples %>%
      #--- need to extract a metric to determine if values are NA ---#
      extract_metrics(mraster = raster[[1]], existing = ., quiet = TRUE) %>%
      #--- remove samples with NA ---#
      stats::na.omit() %>%
      dplyr::select(geometry)

  } else if (location == "centers"){
    
    samples <- sf::st_centroid(gridR) %>%
      dplyr::rename(geometry = x)
    
    #--- check to make sure that samples intersect raster extent (cellsize check) ---#
    if(isFALSE(lengths(suppressMessages(sf::st_intersects(samples, rasterext))) > 0)){
      stop("No samples intersect with 'raster'. Ensure 'cellsize' makes sense.", call. = FALSE)
    }
    
    samples <- samples %>%
      #--- need to extract a metric to determine if values are NA ---#
      extract_metrics(mraster = raster[[1]], existing = ., quiet = TRUE) %>%
      #--- remove samples with NA ---#
      stats::na.omit() %>%
      dplyr::select(geometry)
    
  }
  }
  
  if (isTRUE(plot)) {
    
    #--- plot input raster and random samples ---#
    
    if (!is.null(access)) {
      suppressWarnings(terra::plot(rasterP[[1]]))
      suppressWarnings(terra::plot(gridR, add = TRUE, border = c("blueviolet"), alpha = 0.01))
      suppressWarnings(terra::plot(access_buff$buff, add = T, border = c("gray30"), col = "gray10", alpha = 0.1))
      suppressWarnings(terra::plot(samples, add = TRUE, col = "black"))
    } else {
      suppressWarnings(terra::plot(rasterP[[1]]))
      suppressWarnings(terra::plot(gridR, add = TRUE, border = c("blueviolet"), alpha = 0.01))
      suppressWarnings(terra::plot(samples, add = TRUE, col = "black"))
    }
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
  
  if (isTRUE(details)) {
    
    #--- output metrics details along with stratification raster ---#
    
    output <- list(samples = samples, tessellation = gridR)
    
    #--- output samples dataframe ---#
    return(output)
  } else {
    
    #--- just output raster ---#
    
    return(samples)
  }
}

#--- rotation and translate functions for random rotation ---#
#https://r-spatial.github.io/sf/articles/sf3.html
#https://stackoverflow.com/questions/51282724/creating-a-regular-polygon-grid-over-a-spatial-extent-rotated-by-a-given-angle
rot <-  function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
tran <-  function(geo, ang, center) (geo - center) * rot(ang * pi / 180) + center

