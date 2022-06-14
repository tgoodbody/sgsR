set.seed(2022)
o <- sample_ahels(mraster = mraster, existing = existing, matrices = mat, nSamp = 5, details = TRUE)
o1 <- sample_ahels(mraster = mraster, existing = existing, nQuant = 5, threshold = 0.8, details = TRUE)
o2 <- sample_ahels(mraster = mraster, existing = existing, matrices = mat, threshold = 0.8, tolerance = 0.025)

test_that("Input classes", {
  expect_error(sample_ahels(mraster = "mraster", existing = existing, matrices = mat, nSamp = 5), "'mraster' must be type SpatRaster.")
  expect_error(sample_ahels(mraster = mraster, existing = "existing", matrices = mat, nSamp = 5), "'existing' must be a data.frame or sf object.")
  expect_error(sample_ahels(mraster = mraster, existing = existing.df, matrices = mat, nSamp = "50"), "'nSamp' must be type numeric.")
  expect_error(sample_ahels(mraster = mraster, existing = existing.df, matrices = mat, threshold = "0.9"), "'threshold' must be type numeric.")
  expect_error(sample_ahels(mraster = mraster, existing = existing.df, matrices = mat, threshold = 0.9, tolerance = "0.4"), "'tolerance' must be type numeric.")
  expect_error(sample_ahels(mraster = mraster, existing = existing.df.n, matrices = mat), "'existing' must have columns named 'X' and 'Y'.")
})

test_that("Total outputs", {
  expect_equal(nrow(o$samples), 205L)
  expect_equal(ncol(o$samples), 5L)
  expect_equal(nrow(o$details$existingRatio), 10L)
  expect_equal(nrow(o1$details$existingRatio), 5L)
  expect_equal(nrow(o$samples[o$samples$type == "existing",]), 200L)
  expect_type(o,"list")
  expect_s3_class(o2,"sf")
})

test_that("Messages", {
  expect_message(sample_ahels(mraster = mraster, existing = existing, matrices = mat, threshold = 0.8, tolerance = 0.025), "Threshold of 0.8 with a tolerance of 0.025 provided. Samples will be added until sampling ratios are >= 0.775.")
  expect_message(sample_ahels(mraster = mraster, existing = existing, matrices = mat, nSamp = 5, tolerance = 0.025), "A tolerance of 0.025 has been provided. Samples will be added until 5 is reached or until sampling ratios are all >= 0.975.")
  expect_message(sample_ahels(mraster = mraster, existing = existing), "Creating covariance matrix.")
})