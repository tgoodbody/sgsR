test_that("Total outputs", {
  skip_on_cran()
  
  o <- calculate_coobs(mraster = mraster, existing = existing, plot = TRUE, cores = 4)
  
  expect_equal(nrow(o), 277L)
  expect_equal(ncol(o), 373L)
  expect_equal(nlyr(o), 2L)

  expect_named(o, c("COOB","COOBclass"))
  
  expect_s4_class(o,"SpatRaster")

})

