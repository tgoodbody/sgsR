o <- extract_strata(sraster = sraster, existing = existing)
o1 <- extract_strata(sraster = sraster, existing = existing, data.frame = TRUE)

unq <- unique(o$strata)

xx <- c(sraster, sraster)

test_that("Single breaks classes", {
  expect_error(extract_strata(sraster = "A", existing = existing), "'sraster' must be type SpatRaster.")
  expect_error(extract_strata(sraster = xx, existing = existing), "'sraster' must have a single layer named 'strata'.")
  expect_error(extract_strata(sraster = mraster[[1]], existing = existing), "'sraster' must have a layer named 'strata'.")
  expect_error(extract_strata(sraster = sraster, existing = "existing"), "'existing' must be a data.frame or sf object.")
  expect_error(extract_strata(sraster = sraster, existing = existing, quiet = 1), "'quiet' must be type logical.")
  expect_error(extract_strata(sraster = sraster, existing = access), "'existing' must be an 'sf' object of type 'sfc_POINT' geometry.")
})

test_that("Total outputs", {
  expect_equal(nrow(o), 200L)
  expect_equal(ncol(o), 3L)
  expect_equal(c(1, 2, 3, 4), unq)
  expect_s3_class(o, "sf")

  expect_named(o1, c("X", "Y", "strata", "FID"))
  expect_equal(nrow(o1), 200L)
})


test_that("messages", {
  e <- as.data.frame(sf::st_coordinates(existing))

  names(e) <- c("x", "y")

  expect_message(extract_strata(sraster = sraster, existing = e), "Column coordinate names are lowercase - converting to uppercase.")
  expect_message(extract_strata(sraster = sraster, existing = e), "Column coordinate names are lowercase - converting to uppercase.")

  expect_equal(ncol(extract_strata(sraster = sraster, existing = e)), 2L)
  expect_equal(ncol(extract_strata(sraster = sraster, existing = e, data.frame = TRUE)), 3L)

  expect_message(extract_strata(sraster = sraster, existing = existingna), "16 samples are located where strata values are NA.")
})


test_that("errors", {
  expect_error(extract_strata(sraster = sraster, existing = data.frame(strata = c(1, 2, 3), A = c(1, 2, 3), B = c(1, 2, 3))), "existing' must have columns named 'X' and 'Y'.")
  expect_error(extract_strata(sraster = sraster, existing = data.frame(strata = c(1, 2, 3), X = c(1, 2, 3), Y = c(1, 2, 3))), "'existing' only extracts NA values. Ensure that 'existing' overlaps with 'sraster'.")
})
