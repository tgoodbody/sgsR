o <- mask_access(raster = mraster, access = access, buff_inner = NULL, buff_outer = 100)
o1 <- mask_access(raster = mraster, access = access, buff_inner = 20, buff_outer = 100)

odf <- terra::values(o$rast, dataframe=TRUE)

odfna <- odf[complete.cases(odf),]

test_that("Total outputs", {
  expect_equal(nrow(o$rast), 277L)
  expect_equal(ncol(o$rast), 373L)
  expect_equal(length(odfna), 3L)
  expect_equal(nrow(odfna), 28726L)
  expect_s4_class(o$rast,"SpatRaster")
  expect_s4_class(o$buff,"SpatVector")
})


test_that("messages works", {
  expect_message(mask_access(raster = mraster, access = access, buff_inner = NULL, buff_outer = 100),"An access layer has been provided. An external buffer of 100 m have been applied.")
  expect_message(mask_access(raster = mraster, access = access, buff_inner = 50, buff_outer = 100),"An access layer has been provided. An internal buffer of 50 m and an external buffer of 100 m have been applied.")
})
