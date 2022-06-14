o <- strat_kmeans(mraster = mraster, nStrata = 4, iter = 100, details = TRUE)

odf <- terra::values(o$raster, dataframe=TRUE)

odfna <- odf[complete.cases(odf),]

test_that("Total outputs", {
  expect_message(strat_kmeans(mraster = mraster, nStrata = 4, iter = 100),"K-means being performed on 3 layers with 4 centers.")
  expect_equal(length(o$details$cluster), 91195L)
  expect_equal(length(o$details$centers), 4L)
  expect_equal(ncol(o$raster), 373L)
  expect_equal(length(odfna), 91195L)
  expect_equal(length(unique(odfna)), 4L)
})

test_that("Out classes", {
  expect_s4_class(o$raster,"SpatRaster")
  expect_s3_class(o$details,"kmeans")
  expect_equal(c(1,2,3,4),sort(unique(odfna)))
})