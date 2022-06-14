o <- extract_metrics(mraster = mraster, existing = existing)
o1 <- extract_metrics(mraster = mraster, existing = existing, data.frame = TRUE)

test_that("Total outputs", {
  expect_equal(nrow(o), 200L)
  expect_equal(ncol(o), 5L)
  expect_equal(mean(o$zq90), 14.35025)
  expect_s3_class(o,"sf")
  
  expect_named(o1,c("X","Y","strata","zq90","pzabove2","zsd"))
  expect_equal(nrow(o1),200L)
  expect_equal(ncol(o1),6L)
})



