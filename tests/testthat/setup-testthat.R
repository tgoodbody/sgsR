library(terra)
library(dplyr)

#--- metrics raster ---#
mr <- system.file("extdata", "mraster.tif", package = "sgsR")
mraster <- terra::rast(mr)

#--- metrics raster ---#
mr <- system.file("extdata", "mraster_small.tif", package = "sgsR")
mrastersmall <- terra::rast(mr)

#--- strat raster ---#
sr <- system.file("extdata", "sraster.tif", package = "sgsR")
sraster <- terra::rast(sr)

#--- strat raster stack ---#
sraster2 <- c(sraster, sraster)

#--- existing samples ---#
e <- system.file("extdata", "existing.shp", package = "sgsR")
existing <- sf::st_read(e, quiet = TRUE)

existing_samples <- extract_strata(sraster = sraster, existing = existing) %>%
  extract_metrics(., mraster = mraster)

#--- existing with NA samples ---#
ena <- system.file("extdata", "existingna.shp", package = "sgsR")
existingna <- sf::st_read(ena, quiet = TRUE)

#--- access ---#
a <- system.file("extdata", "access.shp", package = "sgsR")
access <- sf::st_read(a)

#--- read polygon coverage ---#
poly <- system.file("extdata", "inventory_polygons.shp", package = "sgsR")
fri <- sf::st_read(poly)

#--- categorical raster ---#
set.seed(2022)
x <- terra::rast(ncol = terra::ncol(sraster), nrow = terra::nrow(sraster), ext = terra::ext(sraster))

x1 <- terra::rast(ncol = 10, nrow = 10, ext = terra::ext(mraster))

terra::values(x) <- sample(LETTERS[1:5], size = terra::ncell(x), replace = TRUE)
names(x) <- "strata"

crs(x) <- crs(mraster)

xmraster <- c(mraster, x)

#--- logical raster ---#
x2 <- x1 %>%
  terra::setValues(., rep(TRUE, 100))
x2 <- c(x2, x2)
names(x2) <- c("strata1", "strata2")

#--- coordinates ---#
coords <- sf::st_coordinates(existing)

#--- dataframes and NA dataframes ---#
existing.df.n.xy <- existing %>%
  extract_metrics(mraster, .) %>%
  sf::st_drop_geometry(.) %>%
  as.data.frame() %>%
  cbind(., coords)

existing.df.n.xy.lc <- existing %>%
  sf::st_drop_geometry(.) %>%
  as.data.frame() %>%
  cbind(., coords)

names(existing.df.n.xy.lc) <- c("FID", "x", "y")

#--- supply quantile and covariance matrices ---#
mat <- calculate_pop(mraster = mraster)

#-- additional ---#
weights <- c(0.25, 0.25, 0.25, 0.25)

e <- extract_strata(sraster, existing)

existing.df <- data.frame(strata = e$strata)
existing.df.n <- data.frame(name = e$strata)

#--- distance to access ---#

d <- calculate_distance(raster = mraster, access = access)

de <- calculate_distance(raster = mraster, access = access) %>%
  extract_metrics(., existing) %>%
  dplyr::select(-FID)

tmp_file <- file.path(tempdir(), "temp.shp")
tmp_file_df <- file.path(tempdir(), "temp.txt")
tmp_file_rast <- file.path(tempdir(), "temp.tif")
