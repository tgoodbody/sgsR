set.seed(2022)
o <- sample_balanced(mraster = mraster, nSamp = 50)

test_that("Total outputs", {
  expect_equal(nrow(o), 50L)
  expect_equal(ncol(o), 1L)
  expect_s3_class(o,"sf")
})

test_that("Messages", {
  expect_error(sample_clhs(mraster = mraster, nSamp = 320, cost = "A"),"No layer named 'A' exists in 'mraster'.")
})  

