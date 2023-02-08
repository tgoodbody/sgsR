o <- calculate_pcomp(mraster = mraster, nComp = 2)
o1 <- calculate_pcomp(mraster = mraster, nComp = 2, scale = FALSE, center = FALSE)

test_that("inputs", {
  expect_error(calculate_pcomp(mraster = "A", nComp = 4), "'mraster' must be type SpatRaster.")
  expect_error(calculate_pcomp(mraster = mraster, nComp = "4"), "'nComp' must be type numeric.")
  expect_error(calculate_pcomp(mraster = mraster, nComp = 4, center = "TRUE"), "'center' must be type logical.")
  expect_error(calculate_pcomp(mraster = mraster, nComp = 4, scale = "TRUE"), "'scale' must be type logical.")
  expect_error(calculate_pcomp(mraster = mraster, nComp = 4, plot = "TRUE"), "'plot' must be type logical.")
  expect_error(calculate_pcomp(mraster = mraster, nComp = 4), "nComp must be <= the number of layers in 'mraster'.")
})


test_that("Total outputs", {
  expect_equal(nrow(o), 277L)
  expect_equal(ncol(o), 373L)
  expect_equal(max(terra::values(o1), na.rm = TRUE), 14.822163)

  expect_named(o, c("PC1", "PC2"))

  expect_s4_class(o, "SpatRaster")

  expect_type(calculate_pcomp(mraster = mraster, nComp = 2, details = TRUE, plot = TRUE), "list")

  expect_message(calculate_pcomp(mraster = mraster, nComp = 2, filename = file.path(tempdir(), "temp.tif"), overwrite = TRUE), "Output raster written to disc.")
})
