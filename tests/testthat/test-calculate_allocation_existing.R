nSamp <- 48

# Test for correct output when allocation is 'prop'
test_that("calculate_allocation_existing() returns correct output for allocation = 'prop'", {
  output <- calculate_allocation_existing(existing = existing_samples, nSamp = nSamp)
  expect_equal(nrow(output), 4L)
  expect_equal(sum(output$total), 48L)
})

# Test for correct output when allocation is 'optim'
test_that("calculate_allocation_existing() returns correct output for allocation = 'optim'", {
  output <- calculate_allocation_existing(existing = existing_samples, nSamp = nSamp, allocation = "optim", metric = "zq90")
  expect_equal(nrow(output), 4L)
  expect_equal(sum(output$total), 48L)
})

# Test for correct output when allocation is 'equal'
test_that("calculate_allocation_existing() returns correct output for allocation = 'equal'", {
  output <- calculate_allocation_existing(existing = existing_samples, nSamp = nSamp, allocation = "equal")
  expect_equal(nrow(output), 4L)
  expect_equal(sum(output$total), nSamp*4)
})

# Test for correct error message when 'existing' is missing 'strata' layer
test_that("calculate_allocation_existing() throws an error when 'existing' is missing 'strata' layer", {
  existing_wrong <- list(counts = c(20, 30, 50))
  expect_error(calculate_allocation_existing(existing = existing_wrong, nSamp = nSamp),
               "'existing must have a layer named 'strata'")
})

# Test for correct error message when 'nSamp' is not numeric
test_that("calculate_allocation_existing() throws an error when 'nSamp' is not numeric", {
  expect_error(calculate_allocation_existing(existing = existing_samples, nSamp = "50"),
               "'nSamp' must be type numeric")
})

# Test for correct error message when allocation type is not recognized
test_that("calculate_allocation_existing() throws an error when allocation type is not recognized", {
  expect_error(calculate_allocation_existing(existing = existing_samples, nSamp = nSamp, allocation = "unknown"),
               "Unknown allocation type: 'unknown'")
})
