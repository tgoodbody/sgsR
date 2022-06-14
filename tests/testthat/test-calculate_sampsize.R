o <- calculate_sampsize(mraster = mraster, plot = TRUE)
o1 <- calculate_sampsize(mraster = mraster, rse = 0.05)

test_that("Total outputs", {
  expect_equal(nrow(o$nSamp), 123L)
  expect_equal(ncol(o$nSamp), 3L)
  expect_equal(mean(o$nSamp$rse), 0.03)
  expect_s3_class(o$plot,"gg")
  
  expect_equal(sum(o1$nSamp), 188L)
  expect_equal(nrow(o1), 3L)
})
