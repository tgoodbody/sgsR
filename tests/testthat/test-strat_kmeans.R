o <- strat_kmeans(mraster = mraster, nStrata = 4, iter = 100, details = TRUE)

odf <- terra::values(o$raster, dataframe=TRUE)

odfna <- odf[complete.cases(odf),]

test_that("errors", {
  
  expect_error(strat_kmeans(mraster = "mraster", nStrata = 4, iter = 100, details = TRUE),"'mraster' must be type SpatRaster.")
  expect_error(strat_kmeans(mraster = mraster, nStrata = "4", iter = 100, details = TRUE),"'nStrata' must be type numeric.")
  expect_error(strat_kmeans(mraster = mraster, nStrata = 4, iter = "100", details = TRUE),"'iter' must be type numeric.")
  expect_error(strat_kmeans(mraster = mraster, nStrata = 4, iter = 100, algorithm = 4),"'algorithm' must be type character.")
  expect_error(strat_kmeans(mraster = mraster, nStrata = 4, iter = 100, center = "TRUE"),"'center' must be type logical.")
  expect_error(strat_kmeans(mraster = mraster, nStrata = 4, iter = 100, scale = "TRUE"),"'scale' must be type logical.")
  expect_error(strat_kmeans(mraster = mraster, nStrata = 4, iter = 100, plot = "TRUE"),"'plot' must be type logical.")
  expect_error(strat_kmeans(mraster = mraster, nStrata = 4, iter = 100, details = "TRUE"),"'details' must be type logical.")

})

test_that("Total outputs", {
  expect_message(strat_kmeans(mraster = mraster, nStrata = 4, iter = 100, plot = TRUE),"K-means being performed on 3 layers with 4 centers.")
  expect_equal(length(o$details$cluster), 91195L)
  expect_equal(length(o$details$centers), 4L)
  expect_equal(ncol(o$raster), 373L)
  expect_equal(length(odfna), 91195L)
  expect_equal(length(unique(odfna)), 4L)
  
  expect_message(strat_kmeans(mraster = mraster, nStrata = 4, iter = 100, filename = file.path(tempdir(), "temp.tif") , overwrite = TRUE), "Output raster written to disc.")
})

test_that("Out classes", {
  expect_s4_class(o$raster,"SpatRaster")
  expect_s3_class(o$details,"kmeans")
  expect_equal(c(1,2,3,4),sort(unique(odfna)))
})

