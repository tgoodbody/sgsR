test_that("Input classes", {
  expect_error(calculate_sampsize(mraster = "mraster"), "'mraster' must be type SpatRaster.")
  expect_error(calculate_sampsize(mraster = mraster, start = "A"), "'start' must be type numeric.")
  expect_error(calculate_sampsize(mraster = mraster, end = "A"), "'end' must be type numeric.")
  expect_error(calculate_sampsize(mraster = mraster, increment = "A"), "'increment' must be type numeric.")
  expect_error(calculate_sampsize(mraster = mraster, plot = 1), "'plot' must be type logical.")
  expect_error(calculate_sampsize(mraster = x), "'mraster' must contain all numeric values.")
  expect_error(calculate_sampsize(mraster = mraster, rse = "A"), "'rse' must be type numeric.")
  expect_error(calculate_sampsize(mraster = mraster, rse = -1), "'rse' must be > 0.")
  expect_error(calculate_sampsize(mraster = mraster, rse = 0.1), "'rse' must be >= 'start' and <= 'end'.")
  expect_message(calculate_sampsize(mraster = mraster, rse = 0.2, start = 0.1, end = 0.3), "'rse' > 0.15 -  are you sure you want this?")
})

test_that("Total outputs", {
  
  o <- calculate_sampsize(mraster = mraster, plot = TRUE)
  o1 <- calculate_sampsize(mraster = mraster, rse = 0.05)
  
  
  expect_equal(nrow(o$nSamp), 123L)
  expect_equal(ncol(o$nSamp), 3L)
  expect_equal(mean(o$nSamp$rse), 0.03)
  expect_s3_class(calculate_sampsize(mraster = mraster, plot = TRUE)$plot,"gg")
  
  expect_equal(sum(o1$nSamp), 188L)
  expect_equal(nrow(o1), 3L)
  
  expect_message(calculate_sampsize(mraster = mraster, rse = 0.01),regexp = NA)
  
  expect_s3_class(calculate_sampsize(mraster = mraster, rse = 0.01, plot = TRUE)$plot, "gg")
})
