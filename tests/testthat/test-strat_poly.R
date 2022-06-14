#--- load input metrics raster ---#
raster <- system.file("extdata", "sraster.tif", package = "sgsR")
sraster <- terra::rast(raster)

#--- read polygon coverage ---#
poly <- system.file("extdata", "inventory_polygons.shp", package = "sgsR")
fri <- sf::st_read(poly)

#--- stratify polygon coverage ---#
#--- specify polygon attribute to stratify ---#

attribute <- "NUTRIENTS"

#--- specify features within attribute & how they should be grouped ---#
#--- as a single vector ---#

features <- c("poor", "rich", "medium")

o <- strat_poly(
  poly = fri,
  attribute = attribute,
  features = features,
  raster = sraster,
  details = TRUE
)

odf <- terra::values(o$raster, dataframe=TRUE)

odfna <- odf[complete.cases(odf),]

#--- or as multiple lists ---#

g1 <- "poor"
g2 <- c("rich", "medium")

features1 <- list(g1, g2)

o1 <- strat_poly(
  poly = fri,
  attribute = attribute,
  features = features1,
  raster = sraster,
  details = TRUE
)

o1df <- terra::values(o1$raster, dataframe=TRUE)

o1dfna <- o1df[complete.cases(o1df),]


test_that("Total outputs", {
  expect_equal(features, o$lookUp$features)
  expect_equal(ncol(o$raster), 373L)
  expect_equal(nrow(odf), 103321L)
  expect_equal(nrow(o$poly), 524L)
  expect_equal(length(odfna), 89878L)
  expect_equal(length(unique(odfna)), 3L)
  expect_equal(length(unique(o1dfna)), 2L)
})

test_that("Out classes", {
  expect_s4_class(o$poly,"SpatVector")
  expect_s4_class(o$raster,"SpatRaster")
  expect_equal(sort(unique(odfna)),o$lookUp$strata)
  expect_equal(3,o$lookUp$strata[3])
  expect_equal(c(1,2,2),o1$lookUp$strata)
})