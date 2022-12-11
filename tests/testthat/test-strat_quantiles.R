o <- strat_quantiles(mraster = mraster$zq90, nStrata = 4)
o1 <- strat_quantiles(mraster = mraster$zq90,mraster2 = mraster$pzabove2, nStrata = 4, nStrata2 = 4,details = TRUE)

odf <- terra::values(o, dataframe=TRUE)
o1df <- terra::values(o1$raster, dataframe=TRUE)

odfna <- odf[complete.cases(odf),]
o1dfna <- o1df[complete.cases(o1df),]

test_that("Single breaks classes", {
  expect_error(strat_quantiles(mraster = "A", nStrata = 4), "'mraster' must be type SpatRaster.")
  expect_error(strat_quantiles(mraster = mraster, nStrata = "A"), "'nStrata' must be type numeric.")
  expect_error(strat_quantiles(mraster = mraster, nStrata = 4, plot = 2), "'plot' must be type logical.")
  expect_error(strat_quantiles(mraster = mraster, nStrata = 4, samp = "A"), "'samp' must be type numeric.")
  expect_error(strat_quantiles(mraster = mraster, nStrata = 4, details = "TRUE"), "'details' must be type logical.")
  expect_error(strat_quantiles(mraster = mraster, nStrata = 4), "Multiple layers detected in 'mraster'. Please define a singular band to stratify.")
  
  expect_message(strat_quantiles(mraster = mraster$zq90, nStrata = 4, plot = TRUE), regexp = NA)
})

test_that("dual breaks classes", {
  expect_error(strat_quantiles(mraster = "A", nStrata = 4), "'mraster' must be type SpatRaster.")
  expect_error(strat_quantiles(mraster = mraster$zq90, nStrata = 4, mraster2 = "A"), "'mraster2' must be type SpatRaster.")
  
  expect_error(strat_quantiles(mraster = mraster$zq90, nStrata = 4, mraster2 = mraster$pzabove2, nStrata2 = "A"), "'nStrata2' must be type numeric.")
  expect_error(strat_quantiles(mraster = mraster$zq90, nStrata = 4, mraster2 = mraster$zq90, nStrata2 = 4), "'mraster' and 'mraster2' must have different metric names.")
  expect_error(strat_quantiles(mraster = mraster$zq90, nStrata = 4, mraster2 = mraster$pzabove2, nStrata2 = 4, plot = TRUE, samp = 10, details = TRUE), "'samp' must be > 0 <= 1.")
  
  expect_message(strat_quantiles(mraster = mraster$zq90, nStrata = 4, mraster2 = mraster$pzabove2, nStrata2 = 4, plot = TRUE), regexp = NA)
})


test_that("Total outputs", {
  
  expect_message(strat_quantiles(mraster = mraster$zq90, nStrata = 4, filename = file.path(tempdir(), "temp.tif"), overwrite = TRUE), "Output raster written to disc.")
  
  expect_equal(nrow(o), 277L)
  expect_equal(ncol(o), 373L)
  expect_equal(length(odfna), 91195L)
  expect_equal(length(unique(odfna)), 4L)
  expect_equal(length(unique(o1dfna)), 16L)
})

test_that("Out classes", {
  expect_s4_class(o,"SpatRaster")
  expect_s4_class(o1$raster,"SpatRaster")
  expect_equal(sort(unique(o1dfna)),o1$details$class)
  expect_equal(0,o1$details$min_mraster2[1])
  expect_equal(93.900002,o1$details$min_mraster2[16])
  expect_lt(o1$details$min_mraster[4],o1$details$max_mraster[4])
})