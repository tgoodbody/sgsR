library(terra)

#--- metrics raster ---#
mr <- system.file("extdata", "mraster.tif", package = "sgsR")
mraster <- terra::rast(mr)

#--- strat raster ---#
sr <- system.file("extdata", "sraster.tif", package = "sgsR")
sraster <- terra::rast(sr)

#--- existing samples ---#
e <- system.file("extdata", "existing.shp", package = "sgsR")
existing <- sf::st_read(e, quiet = TRUE)

#--- access ---#
a <- system.file("extdata", "access.shp", package = "sgsR")
access <- sf::st_read(a)

#--- read polygon coverage ---#
poly <- system.file("extdata", "inventory_polygons.shp", package = "sgsR")
fri <- sf::st_read(poly)

#--- categorical raster ---#
set.seed(2022)
x <- terra::rast(ncol = terra::ncol(sraster), nrow = terra::nrow(sraster), ext = terra::ext(sraster))

terra::values(x) <- sample(LETTERS[1:5], size = terra::ncell(x), replace = TRUE)
names(x) <- "strata"

#--- coordinates ---#
coords <- sf::st_coordinates(existing)

#--- dataframes and NA dataframes ---#
existing.df <- existing %>% sf::st_drop_geometry(.) %>% as.data.frame() %>% cbind(., coords)
existing.df.n <- data.frame(name = existing.df$X)

#--- supply quantile and covariance matrices ---#
mat <- calculate_pop(mraster = mraster)
