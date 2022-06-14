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
  details = TRUE
)

odf <- terra::values(o$raster, dataframe=TRUE)
odfna <- odf[complete.cases(odf),]

o1df <- terra::values(o1$raster, dataframe=TRUE)
o1dfna <- o1df[complete.cases(o1df),]


test_that("Total outputs", {
  expect_equal(nrow(o$raster), 277L)
  expect_equal(ncol(o$raster), 373L)
  expect_equal(length(odfna), 91195L)
  expect_equal(length(unique(odfna)), 5L)
  expect_equal(length(unique(o1dfna)), 17L)
})

test_that("Out classes", {
  expect_s4_class(o$raster,"SpatRaster")
  expect_s4_class(o1$raster,"SpatRaster")
  expect_s3_class(o$plot,"gg")
  expect_equal(sort(unique(o1dfna)),seq(1,17,1))
  expect_equal(c(3,5,11,18),o$details$breaks)
  expect_equal(c(1,2,5),o1$details$breaks2)
})

