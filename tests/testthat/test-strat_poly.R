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

odf <- terra::values(o$raster, dataframe = TRUE)

odfna <- odf[complete.cases(odf), ]

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

o1df <- terra::values(o1$raster, dataframe = TRUE)

o1dfna <- o1df[complete.cases(o1df), ]

test_that("Input classes", {
  expect_error(strat_poly(poly = "fri", attribute = attribute, features = features, raster = sraster), "'poly' must be an 'sf' object.")
  expect_error(strat_poly(poly = access, attribute = attribute, features = features, raster = sraster), "'poly' geometry type must be 'POLYGON' or 'MULTIPOLYGON'.")
  expect_error(strat_poly(poly = fri, attribute = attribute, features = features, raster = sraster, details = 2), "'details' must be type logical.")
  expect_error(strat_poly(poly = fri, attribute = TRUE, features = features, raster = sraster), "'attribute' must be type character.")
  expect_error(strat_poly(poly = fri, attribute = 2, features = features, raster = sraster), "'attribute' must be type character.")
  expect_error(strat_poly(poly = fri, attribute = attribute, features = 2, raster = sraster), "'attribute' does not have specified 'features'.")
  expect_error(strat_poly(poly = fri, attribute = attribute, features = features, raster = "sraster"), "'raster' must be type SpatRaster.")
  expect_error(strat_poly(poly = fri, attribute = attribute, features = features, raster = sraster, plot = "TRUE"), "'plot' must be type logical.")
  expect_error(strat_poly(poly = fri, attribute = attribute, features = features, raster = sraster, filename = 2), "'filename' must be a file path character string.")
  expect_error(strat_poly(poly = fri, attribute = attribute, features = features, raster = sraster, filename = file.path(tempdir(), "temp.tif"), overwrite = "A"), "'overwrite' must be type logical.")
})


test_that("Total outputs", {
  expect_equal(features, o$lookUp$features)
  expect_equal(ncol(o$raster), 373L)
  expect_equal(nrow(odf), 103321L)
  expect_equal(nrow(o$poly), 524L)
  expect_equal(length(odfna), 89878L)
  expect_equal(length(unique(odfna)), 3L)
  expect_equal(length(unique(o1dfna)), 2L)
  expect_message(strat_poly(poly = fri, attribute = attribute, features = features, raster = sraster, filename = file.path(tempdir(), "temp.tif"), overwrite = TRUE, plot = TRUE), "Output raster written to disc.")
})

test_that("Out classes", {
  expect_s4_class(o$poly, "SpatVector")
  expect_s4_class(o$raster, "SpatRaster")
  expect_equal(sort(unique(odfna)), o$lookUp$strata)
  expect_equal(3, o$lookUp$strata[3])
  expect_equal(c(1, 2, 2), o1$lookUp$strata)

  features2 <- c(NA, "poor", "medium")
  features3 <- c("poor", "poor", "medium")

  expect_message(strat_poly(poly = fri, attribute = attribute, features = features2, raster = sraster), "'features' contains NA. Is this on purpose?")
  expect_error(strat_poly(poly = fri, attribute = attribute, features = features3, raster = sraster), "Repeated within 'features': poor")
  expect_error(strat_poly(poly = fri %>% dplyr::select(-NUTRIENTS), attribute = attribute, features = features2, raster = sraster), "'poly' does not have a layer named NUTRIENTS.")
})
