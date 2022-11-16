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
  sraster = srasterfri,
  sraster2 = sraster
)

o1 <- strat_map(
  sraster = srasterfri,
  sraster2 = sraster,
  stack = TRUE,
  details = TRUE
)

o1df <- terra::values(o1$raster, dataframe=TRUE)

o1dfna <- o1df[complete.cases(o1df),]

test_that("Single breaks classes", {
  expect_error(strat_map(sraster = "A", sraster2 = srasterfri), "'sraster' must be type SpatRaster.")
  expect_error(strat_map(sraster = sraster, sraster2 = "srasterfri"), "'sraster2' must be type SpatRaster.")
  expect_error(strat_map(sraster = sraster, sraster2 = srasterfri, stack = "TRUE"), "'stack' must be type logical.")
  expect_error(strat_map(sraster = sraster, sraster2 = srasterfri, overwrite = "TRUE"), "'overwrite' must be type logical.")
  expect_error(strat_map(sraster = sraster, sraster2 = srasterfri, filename = TRUE), "'filename' must be type character.")
  expect_error(strat_map(sraster = sraster, sraster2 = srasterfri, plot = "TRUE"), "'plot' must be type logical.")
  expect_error(strat_map(sraster = sraster, sraster2 = srasterfri, details = "TRUE"), "'details' must be type logical.")
  expect_error(strat_map(sraster = mraster, sraster2 = srasterfri), "'sraster' must only contain 1 layer. Please subset the layer you would like to use for mapping.")
  expect_error(strat_map(sraster = srasterfri, sraster2 = mraster), "'sraster2' must only contain 1 layer. Please subset the layer you would like to use for mapping.")
  
  expect_error(strat_map(sraster = mraster$zq90, sraster2 = srasterfri), "A layer name containing 'strata' does not exist within 'sraster'.")
  expect_error(strat_map(sraster = srasterfri, sraster2 = mraster$zq90), "A layer name containing 'strata' does not exist within 'sraster2'.")

})

test_that("Total outputs", {
  expect_equal(1, length(names(o)))
  expect_equal(c("strata","strata2","stratamapped"), names(o1$raster))
  expect_equal(ncol(o1$raster), 373L)
  expect_equal(nrow(o1df), 103321L)
  expect_equal(sum(o1$lookUp$strata), 24L)
  expect_equal(length(o1dfna), 3L)
  expect_equal(nrow(o1dfna), 87431L)
  expect_equal("14", o1$lookUp$stratamapped[4])
  
  
  expect_message(strat_map(sraster = sraster, sraster2 = srasterfri, filename = file.path(tempdir(), "temp.tif"), overwrite = TRUE),"Output raster written to disc.")
  expect_message(strat_map(sraster = sraster, sraster2 = srasterfri, stack = TRUE, filename = file.path(tempdir(), "temp.tif"), overwrite = TRUE),"Output stack written to disc.")
})

test_that("Out classes", {
  expect_s4_class(o,"SpatRaster")
  expect_equal(1,o1$lookUp$strata[3])
  expect_equal(c(3,3,1,1,3,3,1,1,2,2,2,2),o1$lookUp$strata)
})

#--- categorical raster ---#
test_that("Categorical", {
  expect_equal(strat_map(sraster = x, sraster2 = x, details = TRUE)$lookUp$stratamapped[1],"DD")
  expect_equal("D2",strat_map(sraster = x, sraster2 = sraster, details = TRUE)$lookUp$stratamapped[1])
  expect_equal("3D",strat_map(sraster = srasterfri, sraster2 = x, details = TRUE)$lookUp$stratamapped[1])
})
