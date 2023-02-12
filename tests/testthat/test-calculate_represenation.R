o <- calculate_representation(sraster = sraster, existing = existing)
o1 <- calculate_representation(sraster = sraster, existing = existing, plot = TRUE)

ex <- extract_strata(sraster, existing)

test_that("Total outputs", {
  expect_error(calculate_representation(sraster = "sraster", existing = existing), "'sraster' must be type SpatRaster.")
  expect_error(calculate_representation(sraster = mraster$pzabove2, existing = existing), "A layer name containing 'strata' does not exist within 'sraster'. Use extract_strata().")
  expect_error(calculate_representation(sraster = sraster, existing = existing, drop = "A"), "'drop' must be type numeric.")
  expect_error(calculate_representation(sraster = sraster, existing = existing, drop = 20), "'drop' must be a numeric value between 0-1.")
  expect_error(calculate_representation(sraster = sraster, existing = existing, drop = -20), "'drop' must be a numeric value between 0-1.")
})


test_that("Total outputs", {
  expect_equal(nrow(o), 4L)
  expect_equal(sum(o$srasterFreq), 1L)
  expect_equal(sum(o$nSamp), 200L)
  expect_equal(nrow(calculate_representation(sraster = sraster, existing = ex)), 4L)
  expect_equal(nrow(calculate_representation(sraster = sraster, existing = existing, drop = 0.2)), 4L)
})
