o <- extract_metrics(mraster = mraster, existing = existing)
o1 <- extract_metrics(mraster = mraster, existing = existing, data.frame = TRUE)

test_that("Total outputs", {
  expect_equal(nrow(o), 200L)
  expect_equal(ncol(o), 5L)
  expect_equal(mean(o$zq90), 14.35025)
  expect_equal(nrow(o1), 200L)
  expect_equal(ncol(o1), 6L)

  expect_s3_class(o, "sf")

  expect_named(o1, c("X", "Y", "zq90", "pzabove2", "zsd", "FID"))
})


test_that("df input", {
  e <- as.data.frame(sf::st_coordinates(existing))

  names(e) <- c("x", "y")

  expect_message(extract_metrics(mraster = mraster, existing = e), "Column coordinate names are lowercase - converting to uppercase.")
  expect_message(extract_metrics(mraster = mraster, existing = e), "Column coordinate names are lowercase - converting to uppercase.")

  expect_equal(ncol(extract_metrics(mraster = mraster, existing = e)), 5L)
  expect_equal(ncol(extract_metrics(mraster = mraster, existing = e, data.frame = TRUE)), 6L)
})


test_that("errors", {
  expect_error(extract_metrics(mraster = "mraster", existing = existing), "'mraster' must be type SpatRaster.")
  expect_error(extract_metrics(mraster = mraster, existing = "existing"), "'existing' must be a data.frame or sf object.")
  expect_error(extract_metrics(mraster = mraster, existing = access), "'existing' must be an 'sf' object of type 'sfc_POINT' geometry.")
  expect_error(extract_metrics(mraster = mraster, existing = data.frame(strata = c(1, 2, 3), A = c(1, 2, 3), B = c(1, 2, 3))), "existing' must have columns named 'X' and 'Y'.")
  expect_error(extract_metrics(mraster = mraster, existing = data.frame(strata = c(1, 2, 3), X = c(1, 2, 3), Y = c(1, 2, 3))), "'existing' only extracts NA values. Ensure that 'existing' overlaps with 'mraster'.")
  expect_error(extract_metrics(mraster = mraster, existing = existing, quiet = 1), "'quiet' must be type logical.")

  expect_error(extract_metrics(mraster = mraster, existing = existing, filename = file.path(tempdir(), "temp.tif"), overwrite = "TRUE", data.frame = TRUE), "'overwrite' must be type logical.")
  expect_error(extract_metrics(mraster = mraster, existing = existing, filename = file.path(tempdir(), "temp.tif"), overwrite = "TRUE", data.frame = FALSE), "'overwrite' must be type logical.")
})


test_that("writes to disc", {
  expect_message(extract_metrics(mraster = mraster, existing = existing, filename = file.path(tempdir(), "temp.shp"), overwrite = TRUE), "Output samples written to disc.")
  expect_message(extract_metrics(mraster = mraster, existing = existing, data.frame = TRUE, filename = file.path(tempdir(), "temp.csv"), overwrite = TRUE), "Output samples written to disc.")
})
