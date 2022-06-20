test_that("multiplication works", {
  expect_error(mask_access(raster = mraster, access = "access", buff_inner = 50, buff_outer = 50),"'access' must be an 'sf' object.")
  expect_error(mask_access(raster = mraster, access = existing, buff_inner = 50, buff_outer = 50),"'access' geometry type must be 'LINESTRING' or 'MULTILINESTRING'.")
  expect_error(mask_access(raster = mraster, access = access, buff_outer = NULL),"'buff_outer' must be provided when 'access' is defined.")
  expect_error(mask_access(raster = mraster, access = access, buff_outer = "3"),"'buff_outer' must be type numeric.")
  
  expect_error(mask_access(raster = mraster, access = access, buff_inner = 50, buff_outer = 3),"'buff_inner' must be < 'buff_outer'.")
})
