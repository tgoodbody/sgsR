# Test the function
test_that("It returns an error when the 'iter' argument is not numeric", {
  expect_error(sample_existing_clhs(existing = existing, nSamp = 3, iter = "a"), "'iter' must be type numeric.")
})

test_that("It returns an error when the 'iter' argument is less than or equal to 0", {
  expect_error(sample_existing_clhs(existing, nSamp = 3, iter = -1), "'iter' must be  >= 0.")
})

test_that("It returns an error when the input raster has less than 2 attributes", {
  expect_error(sample_existing_clhs(existing, nSamp = 3, raster = mraster[[1]]), "At least 2 raster attributes are required")
})

test_that("It returns an error when no attributes are present in 'existing' except for geometry", {
  expect_error(sample_existing_clhs(existing[, "geometry"], nSamp = 3), "At least 2 attributes are required")
})

test_that("It returns an error when the 'cost' argument is not numeric or character", {
  expect_error(sample_existing_clhs(existing, nSamp = 3, cost = TRUE), "'cost' must be either type numeric or character.")
})

test_that("It returns an error when the 'cost' index does not exist within the given object", {
  expect_error(sample_existing_clhs(existing, nSamp = 3, cost = 4), "'cost' index does not exist within 'existing'.")
})

test_that("It returns an error when the given layer name does not exist in the object", {
  expect_error(sample_existing_clhs(existing, nSamp = 3, cost = "layer3"), "No layer named 'layer3' exists")
})

test_that("It returns a list when 'details' is TRUE", {
  expect_type(sample_existing_clhs(existing_samples, nSamp = 3, details = TRUE), "list")
})
