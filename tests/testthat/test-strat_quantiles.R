o <- strat_quantiles(mraster = mraster$zq90, nStrata = 4)
o1 <- strat_quantiles(mraster = mraster$zq90,mraster2 = mraster$pzabove2, nStrata = 4, nStrata2 = 4,details = TRUE)

odf <- terra::values(o, dataframe=TRUE)
o1df <- terra::values(o1$raster, dataframe=TRUE)

odfna <- odf[complete.cases(odf),]
o1dfna <- o1df[complete.cases(o1df),]

test_that("Total outputs", {
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