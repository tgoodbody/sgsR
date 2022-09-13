test_that("mask access", {
  expect_error(mask_access(raster = mraster, access = "access", buff_inner = 50, buff_outer = 50),"'access' must be an 'sf' object.")
  expect_error(mask_access(raster = mraster, access = existing, buff_inner = 50, buff_outer = 50),"'access' geometry type must be 'LINESTRING' or 'MULTILINESTRING'.")
  expect_error(mask_access(raster = mraster, access = access, buff_outer = NULL),"'buff_outer' must be provided when 'access' is defined.")
  expect_error(mask_access(raster = mraster, access = access, buff_outer = "3"),"'buff_outer' must be type numeric.")
  
  expect_error(mask_access(raster = mraster, access = access, buff_inner = 50, buff_outer = 3),"'buff_inner' must be < 'buff_outer'.")
})

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


test_that("mask existing", {
  expect_error(mask_existing(existing = existing, access = "access", buff_inner = 50, buff_outer = 50),"'access' must be an 'sf' object.")
  expect_error(mask_existing(existing = existing, access = existing, buff_inner = 50, buff_outer = 50),"'access' geometry type must be 'LINESTRING' or 'MULTILINESTRING'.")
  expect_error(mask_existing(existing = existing, access = access, buff_outer = NULL),"'buff_outer' must be provided when 'access' is defined.")
  expect_error(mask_existing(existing = existing, access = access, buff_outer = "3"),"'buff_outer' must be type numeric.")
  
  expect_error(mask_existing(existing = existing, access = access, buff_inner = 50, buff_outer = 3),"'buff_inner' must be < 'buff_outer'.")
  
  expect_message(mask_existing(existing = existing, access = access, buff_inner = 50, buff_outer = 200),"Masking resulted in an output of 62 potential sample units.")
  
})

o <- mask_existing(existing = existing, access = access, buff_inner = NULL, buff_outer = 100)
o1 <- mask_existing(existing = existing, access = access, buff_inner = 20, buff_outer = 100)

test_that("Total outputs", {
  expect_equal(nrow(o$samples), 68L)
  expect_equal(ncol(o$samples), 2L)
  expect_equal(nrow(o1$samples), 58L)
  expect_equal(ncol(o1$samples), 2L)
  expect_s3_class(o$samples,"sf")
  expect_s3_class(o1$samples,"sf")
  expect_type(o1,"list")
})


test_that("messages works", {
  expect_message(mask_existing(existing = existing, access = access, buff_inner = NULL, buff_outer = 100),"An access layer has been provided. An external buffer of 100 m have been applied.")
  expect_message(mask_existing(existing = existing, access = access, buff_inner = 50, buff_outer = 100),"Masking resulted in an output of 31 potential sample units.")
})


