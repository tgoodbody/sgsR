#--- create vector breaks ---#
br.zq90 <- c(3, 5, 11, 18)
br.pz2 <- c(20, 40, 60, 80)
br.sd <- c(1, 2, 5)

breaks <- list(br.zq90,br.pz2,br.sd)

o <- strat_breaks(
  mraster = mraster$zq90,
  breaks = br.zq90
)

o1 <- strat_breaks(
  mraster = mraster,
  breaks = breaks,
  map = TRUE,
  details = TRUE,
  plot = TRUE
)

test_that("Single breaks classes", {
  expect_error(strat_breaks(mraster = "A", breaks = br.max), "'mraster' must be type SpatRaster.")
  expect_error(strat_breaks(mraster = mraster, breaks = "A"), "'breaks' must be type numeric.")
  expect_error(strat_breaks(mraster = mraster, breaks = breaks, map = "A"), "'map' must be type logical.")
  expect_error(strat_breaks(mraster = mraster, breaks = breaks, plot = 2), "'plot' must be type logical.")
  expect_error(strat_breaks(mraster = mraster, breaks = breaks, details = 2), "'details' must be type logical.")
  expect_error(strat_breaks(mraster = mraster[[1:2]], breaks = data.frame(a=c(3,5,11,18), b = c(20,40,60,80)), plot = TRUE), "`breaks` must be a list of numeric vectors of the same length as `mraster`.")
  expect_error(strat_breaks(mraster = mraster, breaks = list(br.zq90,br.sd)), "`mraster` and `breaks` must have the same number of layers & objects.")
  expect_error(strat_breaks(mraster = mraster$zq90, breaks = 300), "'breaks' contains values > the maximum corresponding 'mraster' value.")
  expect_error(strat_breaks(mraster = mraster$zq90, breaks = -1), "'breaks' contains values < the minimum corresponding 'mraster' value.")

})

test_that("Total outputs", {
  
  expect_message(strat_breaks(mraster = mraster$zq90, breaks = br.zq90, filename = file.path(tempdir(), "temp.tif"), overwrite = TRUE), "Output raster written to disc.")
  
  expect_equal(nrow(o), 277L)
  expect_equal(ncol(o), 373L)
})

test_that("Out classes", {
  expect_s4_class(o,"SpatRaster")
  expect_s4_class(o1$raster,"SpatRaster")
  expect_s3_class(o1$plot,"gg")
  expect_equal(nlyr(o1$raster),4L)
  expect_equal(c(1,2,3),o1$breaks$strata_zsd[1:3])
})

