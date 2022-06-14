set.seed(2022)
o <- sample_clhs(mraster = mraster, existing = existing, nSamp = 320, details = TRUE)

test_that("Total outputs", {
  expect_equal(nrow(o$samples), 320L)
  expect_equal(ncol(o$samples), 5L)
  expect_equal(length(o$clhs), 320L)
  expect_s3_class(o$samples,"sf")
})

test_that("Messages", {
  expect_error(sample_clhs(mraster = mraster, nSamp = 320, cost = "A"),"No layer named 'A' exists in 'mraster'.")
})  
