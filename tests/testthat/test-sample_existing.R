set.seed(2022)
o <- sample_existing(existing = existing, raster = mraster, nSamp = 50)
o1 <- sample_existing(existing = de, cost = 4, nSamp = 20, plot = TRUE)
o2 <- sample_existing(existing = de, raster = d, cost = 4, nSamp = 20)
o3 <- sample_existing(existing = de, raster = d, cost = 4, nSamp = 20, details = TRUE, plot = TRUE)

test_that("Input classes", {
  expect_error(sample_existing(existing = "existing", raster = mraster, nSamp = 5), "'existing' must be a data.frame or sf object.")
  expect_error(sample_existing(existing = existing, raster = "mraster", nSamp = 5), "'raster' must be type SpatRaster.")
  expect_error(sample_existing(existing = existing, raster = mraster, nSamp = "A"), "'nSamp' must be type numeric.")
  expect_error(sample_existing(existing = existing, raster = mraster, nSamp = 20, iter = "A"), "'iter' must be type numeric.")
  expect_error(sample_existing(existing = existing, raster = mraster, nSamp = 20, iter = -1), "'iter' must be  >= 0.")
  expect_error(sample_existing(existing = existing, raster = mraster, nSamp = 20, details = "TRUE"), "'details' must be type logical.")
  expect_error(sample_existing(existing = existing, raster = mraster, nSamp = 5, cost = "a"), "No layer named 'a' exists in 'raster'.")
  expect_error(sample_existing(existing = existing, raster = mraster, nSamp = 50, plot = "A"), "'plot' must be type logical")
  expect_error(sample_existing(existing = existing, nSamp = 100, cost = TRUE), "'cost' must be either type numeric or character.")
  expect_error(sample_existing(existing = existing, nSamp = 100, cost = 13), "'cost' index doest not exist within 'existing'.")
  expect_error(sample_existing(existing = existing, raster = mraster, nSamp = 100, cost = 13), "'cost' index doest not exist within 'raster'.")
  expect_error(sample_existing(existing = existing, raster = mraster$zq90, nSamp = 100), "At least 2 raster attributes are required to generate a matrix for sub-sampling.")
  expect_error(sample_existing(existing = existing, nSamp = 100), "At least 2 attributes are required to generate a matrix for sub-sampling.")
  expect_error(sample_existing(raster = mraster, nSamp = 100, existing = data.frame(x = c(1,2,3), y = c(1,2,3))), "'nSamp' must be less than the total number of 'existing' samples.")
  expect_error(sample_existing(existing = de, raster = mraster, nSamp = 10, filename = TRUE), "'filename' must be a file path character string.")
  expect_error(sample_existing(existing = de, raster = mraster, nSamp = 10, filename = file.path(tempdir(), "temp.shp"), overwrite = "A"), "'overwrite' must be type logical.")
})

test_that("non sf existing", {
  expect_message(sample_existing(existing = existing, raster = mraster, nSamp = 20, details = TRUE), "'existing' does not contain attributes with the same names as 'raster'. Extracting metrics.")
})

test_that("categorical", {
  expect_s3_class(sample_existing(existing = existing, raster = xmraster, nSamp = 20, details = TRUE, plot = TRUE)$plotcat, "gg")
})


test_that("Access", {
  expect_message(sample_existing(existing = existing, raster = mraster, nSamp = 20, access = access, buff_inner = 50, buff_outer = 200, plot = TRUE), "An access layer has been provided. An internal buffer of 50 m and an external buffer of 200 m have been applied.")
  expect_message(sample_existing(existing = existing, raster = mraster, nSamp = 20, access = access, buff_outer = 200), "An access layer has been provided. An external buffer of 200 m have been applied.")
  expect_message(sample_existing(existing = existing, raster = mraster, nSamp = 20, filename = file.path(tempdir(), "temp.shp"), overwrite = TRUE), "Output samples written to disc.")
  expect_message(sample_existing(raster = mraster, nSamp = 1, existing = existingna), "16 samples are located where metric values are NA.")
  
  expect_message(sample_existing(existing = existing, raster = mraster, nSamp = 50), "Sub-sampling based on 'raster' distributions.")
  expect_message(sample_existing(existing = de, cost = 4, nSamp = 20), "Sub-sampling based on ALL 'existing' metric distributions. Ensure only attributes of interest are included.")
  expect_message(sample_existing(existing = de, cost = 4, nSamp = 20), "Using `dist2access` as sampling constraint.")
})



test_that("Total outputs", {
  expect_equal(nrow(o), 50L)
  expect_equal(ncol(o), 5L)
  expect_s3_class(o,"sf")
  
  expect_equal(nrow(o1), 20L)
  expect_equal(ncol(o1), 5L)
  expect_s3_class(o1,"sf")
  
  expect_equal(nrow(o2), 20L)
  expect_equal(ncol(o2), 5L)
  expect_s3_class(o2,"sf")
  
  expect_type(o3,"list")
  expect_equal(nrow(o3$samples),20L)
  expect_s3_class(o3$plot,"gg")
  expect_s3_class(o3$clhsOut,"cLHS_result")
  
})

