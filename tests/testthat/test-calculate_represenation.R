o <- calculate_representation(sraster = sraster, existing = existing)
o1 <- calculate_representation(sraster = sraster, existing = existing, plot = TRUE)

test_that("Total outputs", {
  expect_equal(nrow(o), 4L)
  expect_equal(sum(o$srasterFreq), 1L)
  expect_equal(sum(o$nSamp), 200L)
  expect_s3_class(o1$plot,"gg")
})
