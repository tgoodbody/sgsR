o <- calculate_pcomp(mraster = mraster, nComp = 2)
o1 <- calculate_pcomp(mraster = mraster, nComp = 2, scale = FALSE)

test_that("Total outputs", {
  expect_equal(nrow(o), 277L)
  expect_equal(ncol(o), 373L)
  expect_named(o, c("PC1","PC2"))
  
  expect_equal(max(terra::values(o1),na.rm=TRUE), 29.279974)
  expect_s4_class(o,"SpatRaster")
})
