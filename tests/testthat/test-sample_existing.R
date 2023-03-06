nSamp <- 40

test_that("sample_existing returns correct number of samples", {
  # test clhs sampling
  samples <- sample_existing(existing_samples, nSamp, type = "clhs")
  expect_equal(nrow(samples), nSamp)

  # test balanced sampling
  samples <- sample_existing(existing_samples, nSamp, type = "balanced")
  expect_equal(nrow(samples), nSamp)

  # test srs sampling
  samples <- sample_existing(existing_samples, nSamp, type = "srs")
  expect_equal(nrow(samples), nSamp)

  # test stratified sampling
  samples <- sample_existing(existing_samples, nSamp, type = "strat")
  expect_equal(nrow(samples), nSamp)
})

test_that("sample_existing raises error for invalid arguments", {
  # test invalid type argument
  expect_error(
    sample_existing(existing_samples, nSamp, type = "invalid"),
    "'type' must be one of 'clhs','balanced', 'srs', 'strat'."
  )
})
