# Test for write_samples function
test_that("filename must be a character string", {
  expect_error(write_samples(data.frame(), 1))
})

test_that("overwrite must be a logical value", {
  expect_error(write_samples(samples = data.frame(), overwrite = "file.txt", filename = "true"))
})

test_that("writing to file works", {
  write_samples(samples = existing, filename = tmp_file, overwrite = TRUE)
  expect_true(file.exists(tmp_file))
})

test_that("existing file with overwrite=FALSE raises an error", {
  expect_error(write_samples(samples = data.frame(), filename = tmp_file, overwrite = FALSE))
})

# Test for write_samples_df function

test_that("filename must be a character string", {
  expect_error(write_samples_df(data.frame(), 1))
})

test_that("overwrite must be a logical value", {
  expect_error(write_samples_df(data.frame(), "file.txt", "true"))
})

test_that("writing to file works", {
  write_samples_df(samples = existing, filename = tmp_file_df, overwrite = TRUE)
  expect_true(file.exists(tmp_file_df))
})

test_that("existing file with overwrite=FALSE raises an error", {
  expect_error(write_samples_df(samples = data.frame(), filename = tmp_file_df, overwrite = FALSE))
})

# Test for write_raster function

test_that("filename must be a character string", {
  expect_error(write_raster(raster = mrastersmall, filename = 1))
})

test_that("overwrite must be a logical value", {
  expect_error(write_raster(raster = mrastersmall, filename = "file.txt", overwrite = "true"))
})

test_that("writing to file works", {
  write_raster(raster = mrastersmall, filename = tmp_file_rast, overwrite = TRUE)
  expect_true(file.exists(tmp_file_df))
})

test_that("existing file with overwrite=FALSE raises an error", {
  expect_error(write_raster(raster = mrastersmall, filename = tmp_file_rast, overwrite = FALSE))
})
