set.seed(2022)
o <- sample_srs(raster = mraster, nSamp = 50)
o1 <- sample_srs(raster = mraster, nSamp = 20, mindist = 200)

test_that("Input classes", {
  expect_error(sample_srs(raster = "mraster", nSamp = 5), "'raster' must be type SpatRaster.")
  expect_error(sample_srs(raster = mraster, nSamp = "A"), "'nSamp' must be type numeric.")
  expect_error(sample_srs(raster = mraster, nSamp = 5, mindist = "A"), "'mindist' must be type numeric.")
  expect_error(sample_srs(raster = mraster, nSamp = 50, plot = "A"), "'plot' must be type logical")
  expect_error(sample_srs(raster = mraster, nSamp = 100, plot = 1), "'plot' must be type logical.")
  expect_error(sample_srs(raster = mraster, nSamp = 320, filename = TRUE), "'filename' must be a file path character string.")
  expect_error(sample_srs(raster = mraster, nSamp = 10, filename = file.path(tempdir(), "temp.shp"), overwrite = "A"), "'overwrite' must be type logical.")
})

test_that("Access", {
  expect_message(sample_srs(raster = mraster, nSamp = 20, access = access, buff_inner = 50, buff_outer = 200, plot = TRUE), "An access layer has been provided. An internal buffer of 50 m and an external buffer of 200 m have been applied.")
  expect_message(sample_srs(raster = mraster, nSamp = 20, access = access, buff_outer = 200), "An access layer has been provided. An external buffer of 200 m have been applied.")
  expect_message(sample_srs(raster = mraster, nSamp = 20),regexp = NA)
  expect_message(sample_srs(raster = mraster, nSamp = 20, filename = file.path(tempdir(), "temp.shp"), overwrite = TRUE), "Output samples written to disc.")
})

test_that("Total outputs", {
  expect_equal(nrow(o), 50L)
  expect_equal(ncol(o), 1L)
  expect_equal(nrow(o1), 20L)
  expect_s3_class(o,"sf")
})

test_that("Total outputs", {
  skip_on_cran()
  set.seed(2023)
  expect_message(sample_srs(raster = mraster, nSamp = 50, mindist = 1000), "Sampling was not able to select 50 sample units. Output has 38 sample units.")
})


