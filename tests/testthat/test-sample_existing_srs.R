# Test if function returns the expected class
test_that("returns a data frame or spatial data frame", {
  s <- sample_existing_srs(existing, 100)
  expect_true(is.data.frame(s) || inherits(s, "sf"))
})

# Test if the function returns the expected number of rows
test_that("returns the expected number of rows", {
  s <- sample_existing_srs(existing, 100)
  expect_equal(nrow(s), 100L)
})


# Test if the function throws an error when filename is not a character string
test_that("throws an error when filename is not a character string", {
  expect_error(sample_existing_srs(existing, 100, filename = 123))
})
