o <- sample_systematic(raster = mraster, cellsize = 1000,square = FALSE, plot = TRUE)
o1 <- sample_systematic(raster = mraster, cellsize = 2000, details = TRUE)

test_that("Total outputs", {
  expect_equal(nrow(o), 47L)
  expect_equal(ncol(o), 1L)
  expect_equal(nrow(o1$samples), 11L)
  expect_s3_class(o,"sf")
  expect_s3_class(o1$tessellation,"sf")
})




