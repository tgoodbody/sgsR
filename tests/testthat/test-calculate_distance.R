o <- calculate_distance(raster = mraster, access)

test_that("inputs", {
  expect_error(calculate_distance(raster = "A", access = access), "'raster' must be type SpatRaster.")
  expect_error(calculate_distance(raster = mraster, access = "access"), "'access' must be an 'sf' object.")
  expect_error(calculate_distance(raster = mraster, access = existing), "'access' geometry type must be 'LINESTRING' or 'MULTILINESTRING'.")

})


test_that("Total outputs", {
  skip_on_cran()

  expect_equal(nrow(o), 277L)
  expect_equal(ncol(o), 373L)
  expect_equal(max(terra::values(o[[4]]),na.rm=TRUE), 1061.65993)
  
  expect_named(o, c("zq90", "pzabove2", "zsd", "dist2access"))
  
  expect_s4_class(o,"SpatRaster")
  
  expect_type(calculate_pcomp(mraster = mraster, nComp = 2, details = TRUE, plot = TRUE), "list")
  
  expect_message(calculate_distance(raster = mraster, access, filename = file.path(tempdir(), "temp.tif"), overwrite = TRUE), "Output raster written to disc.")
})
