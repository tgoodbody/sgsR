o <- extract_strata(sraster = sraster, existing = existing)
o1 <- extract_strata(sraster = sraster, existing = existing, data.frame = TRUE)

unq <- unique(o$strata)

test_that("Total outputs", {
  expect_equal(nrow(o), 200L)
  expect_equal(ncol(o), 2L)
  expect_equal(c(1,2,3,4), unq)
  expect_s3_class(o,"sf")
  
  expect_named(o1,c("X","Y","strata"))
  expect_equal(nrow(o1),200L)
})


