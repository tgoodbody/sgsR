toSample <- calculate_allocation_existing(existing_samples, 100)

test_that("function returns a data frame", {
  result <- sample_existing_strat(existing_samples, toSample, nSamp = 100)
  expect_s3_class(result, "sf")
})

test_that("function returns correct number of rows", {
  result <- sample_existing_strat(existing_samples, toSample, nSamp = 100)
  expected_rows <- sum(toSample$total)
  expect_equal(nrow(result), expected_rows)
})

test_that("function throws error when strata variable not found", {
  expect_error(sample_existing_strat(existing_samples, toSample, nSamp = 100, strata = "invalid"))
})

test_that("function throws error when filename exists and overwrite = FALSE", {
  expect_error(sample_existing_strat(existing_samples, toSample, nSamp = 100, filename = filename, overwrite = FALSE))
})

test_that("function throws error when total sample size exceeds population size in a stratum", {
  invalid_sample <- data.frame(strata = "a", total = 100)
  expect_error(sample_existing_strat(existing_samples, invalid_sample, nSamp = 100))
})

test_that("take_samples returns correct number of samples", {
  samples <- take_samples(existing_samples, toSample, strata = 1)
  expect_equal(nrow(samples), 25L)
})

test_that("take_samples raises error when not enough samples in strata", {
  expect_error(take_samples(existing_samples, toSample = toSample, strata = 199),
               "Not enough samples in strata: 199 to take:  sample units.")

})