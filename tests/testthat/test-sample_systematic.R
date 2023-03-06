test_that("Input classes", {
  skip_on_cran()
  expect_error(sample_systematic(raster = "A", cellsize = 1000), "'raster' must be type SpatRaster.")
  expect_error(sample_systematic(raster = mraster, cellsize = TRUE), "'cellsize' must be type numeric.")
  expect_error(sample_systematic(raster = mraster, cellsize = -1), "'cellsize' must be > 0.")
  expect_error(sample_systematic(raster = mraster, cellsize = 1000, plot = 2), "'plot' must be type logical.")
  expect_error(sample_systematic(raster = mraster, cellsize = 1000, square = "TRUE"), "'square' must be type logical.")
  expect_error(sample_systematic(raster = mraster, cellsize = 1000, location = 3), "'location' must be type character.")
  expect_error(sample_systematic(raster = mraster, cellsize = 1000, location = "not_center"), "'location' must be one of 'centers', 'corners', or 'random'.")
})

test_that("Total outputs", {
  skip_on_cran()
  set.seed(2022)
  o <- sample_systematic(raster = mraster, cellsize = 1000, square = FALSE, plot = TRUE)
  o1 <- sample_systematic(raster = mraster, cellsize = 1000)

  expect_equal(nrow(o), 39L)
  expect_equal(ncol(o), 1L)

  expect_equal(nrow(o1), 40L)
  expect_equal(ncol(o1), 1L)

  expect_equal(nrow(sample_systematic(raster = mraster, cellsize = 2000, details = TRUE)$samples), 7L)
  expect_s3_class(o, "sf")
  expect_s3_class(sample_systematic(raster = mraster, cellsize = 2000, details = TRUE)$tessellation, "sf")
})


test_that("corners", {
  skip_on_cran()
  set.seed(2022)
  o1 <- sample_systematic(raster = mraster, cellsize = 1000, location = "corners")
  o2 <- sample_systematic(raster = mraster, cellsize = 1000, square = FALSE, location = "corners")

  expect_equal(nrow(o1), 195L)
  expect_equal(nrow(o2), 285L)
})

test_that("random", {
  skip_on_cran()
  set.seed(2022)
  or <- sample_systematic(raster = mraster, cellsize = 1000, location = "random")
  or1 <- sample_systematic(raster = mraster, cellsize = 1000, square = FALSE, location = "random")

  expect_equal(nrow(or), 34L)
  expect_equal(nrow(or1), 44L)

  expect_error(sample_systematic(raster = sraster, cellsize = 5000000, location = "random"), "No samples intersect with 'raster'. Ensure 'cellsize' makes sense.")
})


test_that("messages", {
  skip_on_cran()
  set.seed(2022)
  expect_message(sample_systematic(raster = sraster, cellsize = 20000, square = FALSE, location = "random", force = TRUE), "Forcing samples to fall in non NA locations.")
  expect_message(sample_systematic(raster = sraster, cellsize = 5000, location = "random", force = TRUE), "Forcing samples to fall in non NA locations.")
  expect_error(sample_systematic(raster = sraster, cellsize = 5000000, location = "centers"), "No samples intersect with 'raster'. Ensure 'cellsize' makes sense.")
})
