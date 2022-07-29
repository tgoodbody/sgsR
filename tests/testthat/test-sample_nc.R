set.seed(2022)
o <- sample_nc(mraster = mraster, nSamp = 50)
o1 <- sample_nc(mraster = mraster, nSamp = 25, k = 2)

test_that("Input classes", {
  expect_error(sample_nc(mraster = "mraster", nSamp = 5), "'mraster' must be type SpatRaster.")
  expect_error(sample_nc(mraster = mraster, nSamp = "A"), "'nSamp' must be type numeric.")
  expect_error(sample_nc(mraster = mraster, nSamp = 5, k = "A"), "'k' must be type numeric.")
  expect_error(sample_nc(mraster = mraster, nSamp = 4, iter = "100", details = TRUE),"'iter' must be type numeric.")
  expect_error(sample_nc(mraster = mraster, nSamp = 4, iter = 100, algorithm = 4),"'algorithm' must be type character.")
  expect_error(sample_nc(mraster = mraster, nSamp = 50, plot = "A"), "'plot' must be type logical")
  expect_error(sample_nc(mraster = mraster, nSamp = 100, plot = 1), "'plot' must be type logical.")
  expect_error(sample_nc(mraster = mraster, nSamp = 5, filename = TRUE), "'filename' must be a file path character string.")
  expect_error(sample_nc(mraster = mraster, nSamp = 10, filename = file.path(tempdir(), "temp.shp"), overwrite = "A"), "'overwrite' must be type logical.")
})

test_that("Access", {
  expect_message(sample_nc(mraster = mraster, nSamp = 5, access = access, buff_inner = 50, buff_outer = 200, plot = TRUE), "An access layer has been provided. An internal buffer of 50 m and an external buffer of 200 m have been applied.")
  expect_message(sample_nc(mraster = mraster, nSamp = 5, access = access, buff_outer = 200), "An access layer has been provided. An external buffer of 200 m have been applied.")
  expect_message(sample_nc(mraster = mraster, nSamp = 5),"K-means being performed on 3 layers with 5 centers.")
  expect_message(sample_nc(mraster = mraster, nSamp = 5, filename = file.path(tempdir(), "temp.shp"), overwrite = TRUE), "Output samples written to disc.")
})

test_that("Single metric tests", {
  o3 <- sample_nc(mraster = mraster$zq90, nSamp = 20, details = TRUE)
  
  expect_equal(nrow(o3$samples), 20L)
  expect_equal(ncol(o3$samples), 3L)
  
  expect_type(o3, "list")
})

test_that("Total outputs", {
  expect_equal(nrow(o), 50L)
  expect_equal(ncol(o), 5L)
  expect_equal(nrow(o1), 50L)
  expect_s3_class(o,"sf")
})
