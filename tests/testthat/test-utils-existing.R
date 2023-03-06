# tests for check_existing

test_that("X & Y", {
  expect_message(check_existing(existing = existing.df.n.xy.lc, raster = mraster, nSamp = 5), "'existing' column coordinate names are lowercase - converting to uppercase.")
  expect_error(check_existing(existing = existing.df.n, raster = mraster, nSamp = 5), "'existing' must have columns named 'X' and 'Y'.")
})

test_that("existing must be a data.frame or sf object", {
  expect_error(check_existing(1, raster = NULL, nSamp = 1), "'existing' must be a data.frame or sf object.")
  expect_error(check_existing("existing", raster = NULL, nSamp = 1), "'existing' must be a data.frame or sf object.")
})
test_that("raster must be type SpatRaster", {
  expect_error(check_existing(existing, raster = 1, nSamp = 1), "'raster' must be type SpatRaster.")
})
test_that("nSamp must be type numeric of type integer or double and less than number of existing samples", {
  expect_error(check_existing(existing, raster = NULL, nSamp = "a"), "'nSamp' must be type numeric of type integer or double.")
  expect_error(check_existing(existing, raster = NULL, nSamp = 250), "'nSamp' must be less than the total number of 'existing' samples.")
})
test_that("plot must be type logical", {
  expect_error(check_existing(existing, raster = NULL, nSamp = 1, plot = 1), "'plot' must be type logical.")
})
test_that("details must be type logical", {
  expect_error(check_existing(existing, raster = NULL, nSamp = 1, details = 1), "'details' must be type logical.")
})

test_that("prepare_existing returns an sf object when input is not an sf object", {
  expect_true(inherits(prepare_existing(coords_existing(existing)), "data.frame"))
})

test_that("prepare_existing returns input when it is already an sf object", {
  output <- prepare_existing(existing = existing, raster = NULL, access = NULL, buff_inner = NULL, buff_outer = NULL)
  expect_identical(output, existing)
})

test_that("prepare_existing throws an error when input does not have 'X' and 'Y' columns", {
  expect_error(prepare_existing(existing = existing %>% sf::st_drop_geometry(), raster = NULL, access = NULL, buff_inner = NULL, buff_outer = NULL), "'existing' must have columns named 'X' and 'Y'.")
})

test_that("prepare_existing converts lowercase 'x' and 'y' column names to uppercase 'X' and 'Y'", {
  existing_lower_coords <- prepare_existing(coords_existing(existing)) %>% dplyr::rename(x = X, y = Y)
  output <- prepare_existing(existing = existing_lower_coords, raster = NULL, access = NULL, buff_inner = NULL, buff_outer = NULL)
  expect_true(all(c("X", "Y") %in% colnames(output)))
})

test_that("prepare_existing extracts metrics when raster is supplied", {
  output <- prepare_existing(existing = existing, raster = mraster, access = NULL, buff_inner = NULL, buff_outer = NULL)
  expect_true(all(c("zq90", "pzabove2", "zsd", "FID", "geometry") %in% colnames(output)))

  expect_message(prepare_existing(existing = existing, raster = mraster, access = NULL, buff_inner = NULL, buff_outer = NULL), "'existing' does not contain attributes with the same names as 'raster'. Extracting metrics.")

  expect_message(prepare_existing(existing = existing, raster = mraster, access = access, buff_inner = 100, buff_outer = 400), "An access layer has been provided. An internal buffer of 100 m and an external buffer of 400 m have been applied.")
})

test_that("prepare_existing masks samples with access constraint when access object is supplied", {
  expect_error(prepare_existing(existing = existing, raster = mraster, access = access, buff_inner = NULL, buff_outer = NULL), "'buff_outer' must be provided when 'access' is defined.")
  output <- prepare_existing(existing = existing, raster = mraster, access = access, buff_inner = NULL, buff_outer = 300)
  expect_message(prepare_existing(existing = existing, raster = mraster, access = access, buff_inner = NULL, buff_outer = 300), "Masking resulted in an output of 119 potential sample units.")
  expect_equal(nrow(output), 119L)
})

test_that("coords_existing returns a data.frame with X and Y columns", {
  # check if the output of the function is a data.frame
  expect_true(is.data.frame(coords_existing(existing)))

  # check if the output data.frame contains X and Y columns
  expect_true("X" %in% colnames(coords_existing(existing)))
  expect_true("Y" %in% colnames(coords_existing(existing)))
})

test_that("coords_existing preserves all columns from the input sf object", {
  # check if the function preserves all columns from the input sf object
  expect_equal(colnames(coords_existing(existing_samples)), c("X", "Y", colnames(existing_samples)[1:5]))
  expect_equal(nrow(coords_existing(existing_samples)), nrow(existing_samples))
})
