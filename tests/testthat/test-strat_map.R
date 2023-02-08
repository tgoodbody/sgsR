#--- stratify polygon coverage ---#
#--- specify polygon attribute to stratify ---#

attribute <- "NUTRIENTS"

#--- specify features within attribute & how they should be grouped ---#
#--- as a single vector ---#

features <- c("poor", "rich", "medium")

srasterfri <- strat_poly(
  poly = fri,
  attribute = attribute,
  features = features,
  raster = sraster
)

#--- map srasters ---#
o <- strat_map(
  c(srasterfri, sraster)
)

o1 <- strat_map(
  c(srasterfri, sraster),
  stack = TRUE,
  details = TRUE
)

o1df <- terra::values(o1$raster, dataframe = TRUE)

o1dfna <- o1df[complete.cases(o1df), ]

test_that("Single breaks classes", {
  expect_error(strat_map(sraster = "A"), "'sraster' must be type SpatRaster or a list.")
  expect_error(strat_map(sraster = list(sraster)), "List must have at least 2 'SpatRaster' objects.")
  expect_error(strat_map(sraster = sraster, stack = "TRUE"), "'stack' must be type logical.")
  expect_error(strat_map(sraster = sraster, overwrite = "TRUE"), "'overwrite' must be type logical.")
  expect_error(strat_map(sraster = c(sraster, sraster), filename = TRUE), "'filename' must be type character.")
  expect_error(strat_map(sraster = sraster, plot = "TRUE"), "'plot' must be type logical.")
  expect_error(strat_map(sraster = sraster, details = "TRUE"), "'details' must be type logical.")
  expect_error(strat_map(sraster = x1), "'sraster' has no values.")
  expect_error(strat_map(sraster = sraster), "'sraster' must contain at least 2 layers. Please provide a 'SpatRaster' stack or a list of 'SpatRaster' objects.")
  expect_error(strat_map(sraster = mraster), "A layer name containing 'strata' does not exist within 'sraster'.")

  expect_error(strat_map(sraster = c(x2, x2)), "'SpatRaster' layers must be of class 'numeric', 'integer', 'factor', or 'character'.")
})

test_that("Total outputs", {
  expect_equal(1, length(names(o)))
  expect_equal(c("strata_1", "strata_2", "strata"), names(o1$raster))
  expect_equal(ncol(o1$raster), 373L)
  expect_equal(nrow(o1df), 103321L)
  expect_equal(sum(o1$lookUp$strata_1), 24L)
  expect_equal(length(o1dfna), 3L)
  expect_equal(nrow(o1dfna), 87431L)
  expect_equal(14, o1$lookUp$strata[4])

  expect_message(strat_map(sraster = c(srasterfri, sraster), filename = file.path(tempdir(), "temp.tif"), overwrite = TRUE), "Output raster written to disc.")
  expect_message(strat_map(sraster = c(srasterfri, sraster), stack = TRUE, filename = file.path(tempdir(), "temp.tif"), overwrite = TRUE), "Output stack written to disc.")
})

test_that("Out classes", {
  expect_s4_class(o, "SpatRaster")
  expect_equal(1, o1$lookUp$strata_1[3])
  expect_equal(c(rep(1, 4), rep(2, 4), rep(3, 4)), o1$lookUp$strata_1)
})

#--- categorical raster ---#
test_that("Categorical", {
  expect_equal(strat_map(sraster = c(x, x), details = TRUE)$lookUp$strata[1], "AA")
  expect_equal("A1", strat_map(sraster = c(x, sraster), details = TRUE)$lookUp$strata[1])
  expect_equal("1A", strat_map(sraster = c(srasterfri, x), details = TRUE)$lookUp$strata[1])
})
