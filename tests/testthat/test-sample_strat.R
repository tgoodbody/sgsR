o <- sample_strat(sraster = sraster, nSamp = 20)
oo <- sample_strat(sraster = sraster, nSamp = 20, existing = existing, include = TRUE, remove = TRUE)
o1 <- sample_strat(sraster = sraster, nSamp = 20, allocation = "equal", details = TRUE)
o2 <- sample_strat(sraster = sraster, nSamp = 100, allocation = "manual", weights = c(0.2,0.2,0.2,0.4), details = TRUE)
o3 <- sample_strat(sraster = sraster,mraster = mraster$zq90, nSamp = 100, allocation = "optim", weights = c(0.2,0.2,0.2,0.4), details = TRUE, force = TRUE)

test_that("Total outputs", {
  expect_equal(nrow(o), 20L)
  expect_equal(ncol(o), 4L)
  
  expect_equal(nrow(oo), 20L)
  expect_equal(ncol(oo), 4L)
})

test_that("Test equal", {
  expect_message(sample_strat(sraster = sraster,mraster = mraster$zq90, nSamp = 10, allocation = "equal", weights = c(0.2,0.2,0.2,0.4)), "'weights' was specified but 'allocation = equal' - did you mean to use 'allocation = manual'?")
  expect_message(sample_strat(sraster = sraster,mraster = mraster$zq90, nSamp = 10, allocation = "equal", weights = c(0.2,0.2,0.2,0.4)), "'mraster' was specified but 'allocation = equal' - did you mean to use 'allocation = optim'?")
  expect_equal(nrow(o1$sampleDist), 4L)
  expect_equal(sum(o1$sampleDist$total), 80L)
  expect_equal(unique(o1$samples$type), "new")
  
  expect_equal(sum(o2$sampleDist$total), 100L)
  
  expect_type(o2,"list")
  expect_s3_class(o1$samples,"sf")
})

test_that("Test manual", {
  expect_message(sample_strat(sraster = sraster,mraster = mraster$zq90, nSamp = 10, allocation = "manual", weights = c(0.2,0.2,0.2,0.4)), "'mraster' was specified but 'allocation = manual' - did you mean to use 'allocation = optim'?")
  expect_equal(o2$sampleDist$total[1], 20L)
  
  expect_type(o2,"list")
  expect_s3_class(o2$samples,"sf")
})

test_that("Test optim", {
  expect_message(sample_strat(sraster = sraster,mraster = mraster$zq90, nSamp = 10, allocation = "optim", weights = c(0.2,0.2,0.2,0.4), details = TRUE), "'weights' was specified but 'allocation = optim' - did you mean to use 'allocation = manual'?")
  expect_message(sample_strat(sraster = sraster,mraster = mraster$zq90, nSamp = 100, allocation = "optim", weights = c(0.2,0.2,0.2,0.4), details = TRUE),"nSamp of 100 is not perfectly divisible based on strata distribution. nSamp of 99 will be returned. Use 'force = TRUE' to brute force to 100.")
  expect_equal(sum(o3$sampleDist$total), 100L)
  
  expect_equal(sum(o2$sampleDist$total), 100L)
  
  expect_type(o3,"list")
  expect_s3_class(o3$samples,"sf")
})