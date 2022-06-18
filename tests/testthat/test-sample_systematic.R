test_that("Input classes", {
  expect_error(sample_systematic(raster = "A", cellsize = 1000), "'raster' must be type SpatRaster.")
  expect_error(sample_systematic(raster = mraster, cellsize = TRUE), "'cellsize' must be type numeric.")
  expect_error(sample_systematic(raster = mraster, cellsize = -1), "'cellsize' must be > 0.")
  expect_error(sample_systematic(raster = mraster, cellsize = 1000, plot = 2), "'plot' must be type logical.")
  expect_error(sample_systematic(raster = mraster, cellsize = 1000, square = "TRUE"), "'square' must be type logical.")
  expect_error(sample_systematic(raster = mraster, cellsize = 1000, location = 3), "'location' must be type character.")
  expect_error(sample_systematic(raster = mraster, cellsize = 1000, location = "not_center"), "'location' must be one of 'centers', 'corners', or 'random'.")
  expect_error(sample_systematic(raster = mraster, cellsize = 1000, filename = file.path(tempdir(), "temp.shp"), overwrite = "A"), "'overwrite' must be type logical.")
  expect_error(sample_systematic(raster = mraster, cellsize = 1000, filename = 3, overwrite = "A"), "'filename' must be a file path character string.")
})

test_that("Total outputs", {
  o <- sample_systematic(raster = mraster, cellsize = 1000,square = FALSE, plot = TRUE)
  o1 <- sample_systematic(raster = mraster, cellsize = 1000)
  
  expect_equal(nrow(o), 47L)
  expect_equal(ncol(o), 1L)
  
  expect_equal(nrow(o1), 40L)
  expect_equal(ncol(o1), 1L)
  
  expect_equal(nrow(sample_systematic(raster = mraster, cellsize = 2000, details = TRUE)$samples), 11L)
  expect_s3_class(o,"sf")
  expect_s3_class(sample_systematic(raster = mraster, cellsize = 2000, details = TRUE)$tessellation,"sf")
})


test_that("corners", {
  expect_equal(nrow(sample_systematic(raster = mraster, cellsize = 1000, location = "corners")), 45L)
  expect_equal(nrow(sample_systematic(raster = mraster, cellsize = 1000, square = FALSE, location = "corners")), 312L)
})

test_that("random", {
  set.seed(2022)
  or <- sample_systematic(raster = mraster, cellsize = 1000, location = "random")
  or1 <- sample_systematic(raster = mraster, cellsize = 1000, square = FALSE, location = "random")
  
  expect_equal(nrow(or), 37L)
  expect_equal(nrow(or1), 43L)
})


test_that("messages", {
  expect_message(sample_systematic(raster = sraster, cellsize = 2000, access = access, buff_inner = 50, buff_outer = 200, plot = TRUE), "An access layer has been provided. An internal buffer of 50 m and an external buffer of 200 m have been applied.")
  expect_message(sample_systematic(raster = sraster, cellsize = 2000, access = access, buff_outer = 200), "An access layer has been provided. An external buffer of 200 m have been applied.")
  expect_message(sample_systematic(raster = sraster, cellsize = 2000, access = access, buff_outer = 200, filename = file.path(tempdir(), "temp.shp") , overwrite = TRUE), "Output samples written to disc.")
})

test_that("messages", {
  expect_message(sample_systematic(raster = sraster, cellsize = 200000, square = FALSE, location = "random", force = TRUE),"Forcing samples to fall in non NA locations.")
  expect_message(sample_systematic(raster = sraster, cellsize = 200000, location = "random", force = TRUE),"Forcing samples to fall in non NA locations.")
})
