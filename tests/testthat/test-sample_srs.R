set.seed(2022)
o <- sample_srs(raster = mraster, nSamp = 50)
o1 <- sample_srs(raster = mraster, nSamp = 20, mindist = 200)

test_that("Total outputs", {
  expect_equal(nrow(o), 50L)
  expect_equal(ncol(o), 1L)
  expect_equal(nrow(o1), 20L)
  expect_s3_class(o,"sf")
})
