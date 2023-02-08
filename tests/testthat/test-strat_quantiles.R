o <- strat_quantiles(mraster = mraster$zq90, nStrata = 4)
o1 <- strat_quantiles(mraster = mraster, nStrata = list(c(0.1, 0.2), 4, 2), details = TRUE, plot = TRUE, map = TRUE)

test_that("Single breaks classes", {
  expect_error(strat_quantiles(mraster = "A", nStrata = 4), "'mraster' must be type SpatRaster.")
  expect_error(strat_quantiles(mraster = mraster[[1]], nStrata = "A"), "'nStrata' must be type numeric.")
  expect_error(strat_quantiles(mraster = mraster, nStrata = 4, plot = 2), "'plot' must be type logical.")
  expect_error(strat_quantiles(mraster = mraster, nStrata = 4, map = "A"), "'map' must be type logical.")
  expect_error(strat_quantiles(mraster = mraster, nStrata = 4, details = "TRUE"), "'details' must be type logical.")
  expect_error(strat_quantiles(mraster = mraster, nStrata = 4), "`mraster` and `nStrata` must have the same number of layers & objects.")
  expect_error(strat_quantiles(mraster = mraster[[1:2]], nStrata = data.frame(a = c(3, 5, 11, 18), b = c(20, 40, 60, 80))), "`nStrata` must be a list of numeric vectors of the same length as `mraster`.")

  expect_message(strat_quantiles(mraster = mraster$zq90, nStrata = 4, plot = TRUE), regexp = NA)
})

test_that("Total outputs", {
  expect_message(strat_quantiles(mraster = mraster$zq90, nStrata = 4, filename = file.path(tempdir(), "temp.tif"), overwrite = TRUE), "Output raster written to disc.")

  expect_equal(nrow(o), 277L)
  expect_equal(ncol(o), 373L)
})

test_that("Out classes", {
  expect_s4_class(o, "SpatRaster")
  expect_s4_class(o1$raster, "SpatRaster")
  expect_equal(3L, o1$details$strata[3])
  expect_equal(4L, terra::nlyr(o1$raster))
  expect_s3_class(o1$plot, "gg")
})
