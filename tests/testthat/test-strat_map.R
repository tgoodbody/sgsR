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

test_that("Total outputs", {
  expect_equal(1, length(names(o)))
  expect_equal(c("strata","strata2","stratamapped"), names(o1$raster))
  expect_equal(ncol(o1$raster), 373L)
  expect_equal(nrow(o1df), 103321L)
  expect_equal(sum(o1$lookUp$strata), 24L)
  expect_equal(length(o1dfna), 3L)
  expect_equal(nrow(o1dfna), 87431L)
  expect_equal("14", o1$lookUp$stratamapped[4])
})

test_that("Out classes", {
  expect_s4_class(o,"SpatRaster")
  expect_equal(1,o1$lookUp$strata[3])
  expect_equal(c(3,3,1,1,3,3,1,1,2,2,2,2),o1$lookUp$strata)
})

#--- categorical raster ---#
o2 <- strat_map(
  sraster = srasterfri,
  sraster2 = x,
  details = TRUE
)

test_that("Categorical", {
  expect_message(strat_map(sraster = srasterfri, sraster2 = x),"'sraster2' has factor values. Converting to allow mapping.")
  expect_message(strat_map(sraster = x, sraster2 = srasterfri),"'sraster' has factor values. Converting to allow mapping.")
  expect_equal("D_3",strat_map(sraster = x, sraster2 = srasterfri, details = TRUE)$lookUp$stratamapped_cat[1])
  expect_equal("3_D",o2$lookUp$stratamapped_cat[1])
})
