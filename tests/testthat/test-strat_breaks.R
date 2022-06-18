#--- create vector breaks ---#
br.max <- c(3, 5, 11, 18)
br.sd <- c(1, 2, 5)

o <- strat_breaks(
  mraster = mraster$zq90,
  breaks = br.max,
  plot = TRUE,
  details = TRUE
)

o1 <- strat_breaks(
  mraster = mraster$zq90,
  mraster2 = mraster$zsd,
  breaks = br.max,
  breaks2 = br.sd,
  details = TRUE,
  plot = TRUE
)

odf <- terra::values(o$raster, dataframe=TRUE)
odfna <- odf[complete.cases(odf),]

o1df <- terra::values(o1$raster, dataframe=TRUE)
o1dfna <- o1df[complete.cases(o1df),]

test_that("Single breaks classes", {
  expect_error(strat_breaks(mraster = "A", breaks = br.max), "'mraster' must be type SpatRaster.")
  expect_error(strat_breaks(mraster = mraster, breaks = "A"), "'breaks' must be type numeric.")
  expect_error(strat_breaks(mraster = mraster, breaks = br.max, plot = 2), "'plot' must be type logical.")
  expect_error(strat_breaks(mraster = mraster, breaks = br.max, details = 2), "'details' must be type logical.")
  expect_error(strat_breaks(mraster = mraster, breaks = br.max), "Multiple layers detected in 'mraster'. Please define a singular band to stratify.")
  
  expect_error(strat_breaks(mraster = mraster$zq90, breaks = 300), "'breaks' contains values > the maximum 'mraster' value.")

})

test_that("Dual breaks classes", {
  expect_error(strat_breaks(mraster = mraster$zq90, breaks = br.max, mraster2 = "A"), "'mraster2' must be type SpatRaster.")
  expect_error(strat_breaks(mraster = mraster$zq90, breaks = br.max, mraster2 = mraster$pzabove2, breaks2 = "breaks2"), "'breaks2' must be type numeric.")
  expect_error(strat_breaks(mraster = mraster$zq90, breaks = br.max, mraster2 = mraster, breaks2 = "breaks2"), "Multiple layers detected in 'mraster2'. Please define a singular band to stratify.")
  
  expect_error(strat_breaks(mraster = mraster$zq90, breaks = br.max, mraster2 = mraster$pzabove2, breaks2 = 400), "'breaks2' contains values > the maximum 'mraster2' value.")
  
})


test_that("Total outputs", {
  
  expect_message(strat_breaks(mraster = mraster$zq90, breaks = br.max, filename = file.path(tempdir(), "temp.tif"), overwrite = TRUE), "Output raster written to disc.")
  
  expect_equal(nrow(o$raster), 277L)
  expect_equal(ncol(o$raster), 373L)
  expect_equal(length(odfna), 91195L)
  expect_equal(length(unique(odfna)), 5L)
  expect_equal(length(unique(o1dfna)), 17L)
})

test_that("Out classes", {
  expect_message(strat_breaks(mraster = mraster$zq90, breaks = br.max),regexp = NA)
  expect_s4_class(o$raster,"SpatRaster")
  expect_s4_class(o1$raster,"SpatRaster")
  expect_s3_class(o$plot,"gg")
  expect_s3_class(o1$plot,"gg")
  expect_equal(sort(unique(o1dfna)),seq(1,17,1))
  expect_equal(c(3,5,11,18),o$details$breaks)
  expect_equal(c(1,2,5),o1$details$breaks2)
})

