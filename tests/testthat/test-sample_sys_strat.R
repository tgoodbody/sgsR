test_that("sample_sys_strat function input validation", {
  # Check error message when sraster is not of type SpatRaster
  expect_error(sample_sys_strat(sraster = "sraster", cellsize = 1), "'sraster' must be type SpatRaster.")

  # Check error message when cellsize is not of type numeric
  expect_error(sample_sys_strat(sraster = sraster, cellsize = "a"), "'cellsize' must be type numeric.")

  # Check error message when cellsize is less than 0
  expect_error(sample_sys_strat(sraster = sraster, cellsize = -1), "'cellsize' must be > 0.")

  # Check error message when plot is not of type logical
  expect_error(sample_sys_strat(sraster = sraster, cellsize = 1, plot = 1), "'plot' must be type logical.")

  # Check error message when square is not of type logical
  expect_error(sample_sys_strat(sraster = sraster, cellsize = 1, square = "a"), "'square' must be type logical.")

  # Check error message when details is not of type logical
  expect_error(sample_sys_strat(sraster = sraster, cellsize = 1, details = "a"), "'details' must be type logical.")

  # Check error message when location is not of type character
  expect_error(sample_sys_strat(sraster = sraster, cellsize = 1, location = 1), "'location' must be type character.")

  expect_error(sample_sys_strat(sraster = sraster, cellsize = 1000, location = "not_center"), "'location' must be one of 'centers', 'corners', or 'random'.")
  expect_error(sample_sys_strat(sraster = sraster, cellsize = 1000, filename = file.path(tempdir(), "temp.shp"), overwrite = "A"), "'overwrite' must be type logical.")
  expect_error(sample_sys_strat(sraster = sraster, cellsize = 1000, filename = 3, overwrite = "A"), "'filename' must be a file path character string.")


  s <- sraster
  names(s) <- "a"

  expect_error(sample_sys_strat(sraster = s, cellsize = 1), "'sraster' must have a layer named `strata`.")
})


test_that("Total outputs", {
  skip_on_cran()
  set.seed(2022)
  o <- sample_sys_strat(sraster = sraster, cellsize = 1000, square = FALSE, plot = TRUE)
  o1 <- sample_sys_strat(sraster = sraster, cellsize = 1000, details = TRUE)

  expect_equal(nrow(o), 40L)
  expect_equal(ncol(o), 2L)

  expect_equal(nrow(o1$samples), 41L)
  expect_equal(ncol(o1$samples), 2L)

  expect_s3_class(o, "sf")
  expect_s3_class(o1$tessellation[[1]], "sf")

  expect_equal(37L, nrow(sample_sys_strat(sraster = x, cellsize = 1000, plot = TRUE)))
})
