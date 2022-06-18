set.seed(2022)
o <- sample_ahels(mraster = mraster, existing = existing, matrices = mat, nSamp = 5, details = TRUE)
o1 <- sample_ahels(mraster = mraster, existing = existing, nQuant = 5, threshold = 0.8, details = TRUE)
o2 <- sample_ahels(mraster = mraster, existing = existing, matrices = mat, threshold = 0.8, tolerance = 0.025)


test_that("Input classes", {
  expect_error(sample_ahels(mraster = "mraster", existing = existing, matrices = mat, nSamp = 5), "'mraster' must be type SpatRaster.")
  expect_error(sample_ahels(mraster = mraster, existing = "existing", matrices = mat, nSamp = 5), "'existing' must be a data.frame or sf object.")
  expect_error(sample_ahels(mraster = mraster, existing = existing, nQuant = "10", nSamp = 5), "'nQuant' must be type numeric.")
  expect_error(sample_ahels(mraster = mraster, existing = existing, threshold = 10, nSamp = 5), "'threshold' must be > 0 and < 1.")
  expect_error(sample_ahels(mraster = mraster, existing = existing, nSamp = 5, details = "TRUE"), "'details' must be type logical.")
  expect_error(sample_ahels(mraster = mraster, existing = existing, matrices = mat, nSamp = "50"), "'nSamp' must be type numeric.")
  expect_error(sample_ahels(mraster = mraster, existing = existing, tolerance = 0.5, threshold = 0.4), "'tolerance' cannot be >= `threshold`.")
  expect_error(sample_ahels(mraster = mraster, existing = existing, matrices = mat, threshold = "0.9"), "'threshold' must be type numeric.")
  expect_error(sample_ahels(mraster = mraster, existing = existing, tolerance = 0.2), "'tolerance' must be > 0 and <= 0.1.")
  expect_error(sample_ahels(mraster = mraster, existing = existing, matrices = 4), "'matrices' must be the output from 'calculate_pop\\(\\)\\'.")
  expect_error(sample_ahels(mraster = mraster, existing = existing, matrices = mat, nQuant = 12), "Number of quantiles in provided 'matrices' does not match 'nQuant'.")
  expect_error(sample_ahels(mraster = sraster, existing = existing, matrices = mat), "'mraster' used to generate 'matrices' must be identical.")
  expect_error(sample_ahels(mraster = mraster, existing = existing.df.n, matrices = mat), "'existing' must have columns named 'X' and 'Y'.")
  
  expect_error(sample_ahels(mraster = mraster, matrices = mat, existing = data.frame(x = c(1,2,3), y = c(1,2,3))), "'existing' only extracts NA values. Ensure that 'existing' overlaps with 'mraster'.")
  expect_error(sample_ahels(mraster = mraster, matrices = mat, existing = existing, filename = 56), "'filename' must be a file path character string.")
  expect_error(sample_ahels(mraster = mraster, matrices = mat, existing = existing, filename = file.path(tempdir(), "temp.shp"), overwrite = "A"), "'overwrite' must be type logical.")
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
  expect_message(sample_ahels(mraster = mraster, matrices = mat, existing = existing.df.n.xy.lc), "'existing' column coordinate names are lowercase - converting to uppercase.")
  expect_message(sample_ahels(mraster = mraster, matrices = mat, existing = existingna), "16 samples are located where metric values are NA.")
  
  expect_message(sample_ahels(mraster = mraster, existing = existing, matrices = mat, threshold = 0.8, tolerance = 0.025), "Threshold of 0.8 with a tolerance of 0.025 provided. Samples will be added until sampling ratios are >= 0.775.")
  expect_message(sample_ahels(mraster = mraster, existing = existing, matrices = mat, nSamp = 5, tolerance = 0.025), "A tolerance of 0.025 has been provided. Samples will be added until 5 is reached or until sampling ratios are all >= 0.975.")
  expect_message(sample_ahels(mraster = mraster[[1]], existing = existing), "Creating covariance matrix.")
  expect_message(sample_ahels(mraster = mraster, matrices = mat, existing = existingna, filename = file.path(tempdir(), "temp.shp") , overwrite = TRUE), "Output samples written to disc.")
})
