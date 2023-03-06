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

  expect_error(calculate_allocation_existing(existing = existing_samples, nSamp = nSamp, allocation = "optim"), "'metric' must be supplied if 'allocation = optim'.")
  expect_error(calculate_allocation_existing(existing = existing_samples, nSamp = nSamp, allocation = "optim", metric = c("zq90", "zsd")), "Multiple character strings detected in 'metric'. Please define a singular metric for allocation.")
  expect_error(calculate_allocation_existing(existing = existing_samples, nSamp = nSamp, allocation = "optim", metric = c("layer1")), "No column named layer1 in 'existing'.")
})

# Test for correct output when allocation is 'equal'
test_that("calculate_allocation_existing() returns correct output for allocation = 'equal'", {
  output <- calculate_allocation_existing(existing = existing_samples, nSamp = nSamp, allocation = "equal")
  expect_equal(nrow(output), 4L)
  expect_equal(sum(output$total), nSamp * 4)
})

# Test for correct error message when 'existing' is missing 'strata' layer
test_that("calculate_allocation_existing() throws an error when 'existing' is missing 'strata' layer", {
  existing_wrong <- list(counts = c(20, 30, 50))
  expect_error(
    calculate_allocation_existing(existing = existing_wrong, nSamp = nSamp),
    "'existing must have a layer named 'strata'"
  )
})

# Test for correct error message when 'nSamp' is not numeric
test_that("calculate_allocation_existing() throws an error when 'nSamp' is not numeric", {
  expect_error(
    calculate_allocation_existing(existing = existing_samples, nSamp = "50"),
    "'nSamp' must be type numeric"
  )
})

# Test for correct error message when allocation type is not recognized
test_that("calculate_allocation_existing() throws an error when allocation type is not recognized", {
  expect_error(
    calculate_allocation_existing(existing = existing_samples, nSamp = nSamp, allocation = "unknown"),
    "Unknown allocation type: 'unknown'"
  )
})

# Test for correct error message when allocation type is not recognized
test_that("check force", {
  expect_error(
    calculate_allocation_existing(existing = existing_samples, nSamp = nSamp, force = "1"),
    "'force' must be type logical."
  )
})



# test manual

# Define test cases
test_that("allocate_existing_manual() raises error if 'weights' is not defined", {
  expect_error(allocate_existing_manual(existing_samples, nSamp = nSamp, weights = NULL, allocation = "manual"))
})

test_that("allocate_existing_manual() raises error if 'weights' is not numeric", {
  expect_error(allocate_existing_manual(existing_samples, nSamp = nSamp, allocation = "manual", weights = "not_numeric"))
})

test_that("allocate_existing_manual() raises error if 'weights' does not add up to 1", {
  expect_error(allocate_existing_manual(existing_samples, nSamp = nSamp, allocation = "manual", weights = c(0.5, 0.6)))
})

test_that("allocate_existing_manual() returns correct allocation weights", {
  expected_output <- data.frame(strata = c(1, 2, 3, 4), total = rep(12, 4))
  output <- allocate_existing_manual(existing_samples, nSamp = nSamp, weights = rep(0.25, 4))
  expect_equal(output$total, expected_output$total)
})


test_that("weights messages", {
  expect_message(calculate_allocation_existing(existing_samples, nSamp = nSamp, allocation = "optim", metric = "zq90", weights = rep(0.25, 4)), "'weights' was specified but 'allocation = optim' - did you mean to use 'allocation = manual'?")
  expect_message(calculate_allocation_existing(existing_samples, nSamp = nSamp, allocation = "equal", weights = rep(0.25, 4)), "'weights' was specified but 'allocation = equal' - did you mean to use 'allocation = manual'?")
  expect_message(calculate_allocation_existing(existing_samples, nSamp = nSamp, allocation = "equal", force = TRUE), "`force = TRUE` has no effect when `allocation = equal'. Ignoring.")
  expect_message(calculate_allocation_existing(existing_samples, nSamp = nSamp, weights = rep(0.25, 4)), "'weights' was specified but 'allocation = prop' - did you mean to use 'allocation = manual'?")

  expect_message(calculate_allocation_existing(existing_samples, nSamp = 101), "nSamp of 101 is not perfectly divisible based on strata distribution. nSamp of 100 will be returned. Use 'force = TRUE' to brute force to 101.")
  expect_message(calculate_allocation_existing(existing_samples, nSamp = 101, force = TRUE), "Forcing 101 total samples.")
})
