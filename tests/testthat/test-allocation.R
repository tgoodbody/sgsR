test_that("input classes", {
  expect_error(calculate_allocation(sraster = "error", nSamp = 100), "'sraster' must be type SpatRaster.")
  expect_error(calculate_allocation(sraster = sraster, nSamp = 100, existing = "error"), "'existing' must be a data.frame or sf object.")
  expect_error(calculate_allocation(sraster = sraster, nSamp = 100, existing = existing.df.n), "'existing' must have an attribute named 'strata'. Consider using extract_strata().")
})

test_that("output df names", {
  expect_named(calculate_allocation(sraster = sraster, nSamp = 100, existing = existing.df), c("strata", "total", "need"))
  expect_named(calculate_allocation(sraster = sraster, nSamp = 100), c("strata", "total"))
})

test_that("Expect output", {
  expect_true(nrow(calculate_allocation(sraster = sraster, nSamp = 100)) == 4)
})


test_that("proportional allocation", {
  expect_equal(sum(calculate_allocation(sraster = sraster, nSamp = 100)$total), 100L)

  expect_message(calculate_allocation(sraster = sraster, nSamp = 100), "Implementing proportional allocation of samples.")

  expect_message(calculate_allocation(sraster = sraster, nSamp = 101, force = TRUE), "Forcing 101 total samples.")
})

test_that("optimal allocation", {
  expect_equal(sum(calculate_allocation(sraster = sraster, allocation = "optim", mraster = mraster$zq90, nSamp = 100, force = TRUE)$total), 100L)

  expect_message(calculate_allocation(sraster = sraster, allocation = "optim", mraster = mraster$zq90, nSamp = 100, force = TRUE), "Implementing optimal allocation of samples based on variability of 'zq90'")

  expect_error(calculate_allocation(sraster = sraster, allocation = "optim", nSamp = 100), "'mraster' must be supplied if 'allocation = optim'.")
  expect_error(calculate_allocation(sraster = sraster, allocation = "optim", mraster = mraster, nSamp = 100), "Multiple layers detected in 'mraster'. Please define a singular band for allocation.")
})

test_that("manual allocation", {
  expect_equal(sum(calculate_allocation(sraster = sraster, allocation = "manual", weights = weights, nSamp = 100)$total), 100L)

  expect_message(calculate_allocation(sraster = sraster, allocation = "manual", weights = weights, nSamp = 100), "Implementing allocation of samples based on user-defined weights.")

  expect_error(calculate_allocation(sraster = sraster, allocation = "manual", weights = 1, nSamp = 100), "'weights' must be the same length as the number of strata in 'sraster'.")
  expect_error(calculate_allocation(sraster = sraster, allocation = "manual", weights = c(1, 2, 3, 4), nSamp = 100), "'weights' must add up to 1.")
})

test_that("equal allocation", {
  expect_equal(sum(calculate_allocation(sraster = sraster, allocation = "equal", nSamp = 100)$total), 400L)

  expect_message(calculate_allocation(sraster = sraster, allocation = "equal", nSamp = 100), "Implementing equal allocation of samples.")
  expect_message(calculate_allocation(sraster = sraster, allocation = "equal", nSamp = 100, force = TRUE), "`force = TRUE` has no effect when `allocation = equal'. Ignoring.")
})

test_that("Check error for incorrect class inputs", {
  expect_error(calculate_allocation(sraster = sraster, nSamp = "A"))
  expect_error(calculate_allocation(sraster = mraster, nSamp = 100))
  expect_error(calculate_allocation(sraster = sraster, nSamp = 100, allocation = "error"))
  expect_error(calculate_allocation(sraster = sraster, nSamp = 100, allocation = "manual", weights = "A"))
  expect_error(calculate_allocation(sraster = sraster, nSamp = 100, allocation = "manual", weights = c(1, 2, 3, 4)))
  expect_error(calculate_allocation(sraster = sraster, nSamp = 100, existing = TRUE))
  expect_error(calculate_allocation(sraster = sraster, nSamp = 100, force = "error"))
})

test_that("Check error for incorrect class inputs", {
  expect_message(calculate_allocation(sraster = sraster, nSamp = 100, weights = weights), "'weights' was specified but 'allocation = prop' - did you mean to use 'allocation = manual'?")
  expect_message(calculate_allocation(sraster = sraster, nSamp = 100, mraster = mraster), "'mraster' was specified but 'allocation = prop' - did you mean to use 'allocation = optim'?")
})
