o <- sample_strat(sraster = sraster, nSamp = 20)
o1 <- sample_strat(sraster = sraster, nSamp = 20, allocation = "equal", details = TRUE)
o2 <- sample_strat(sraster = sraster, nSamp = 100, allocation = "manual", weights = c(0.2,0.2,0.2,0.4), details = TRUE)
o3 <- sample_strat(sraster = sraster,mraster = mraster$zq90, nSamp = 100, allocation = "optim", weights = c(0.2,0.2,0.2,0.4), details = TRUE, force = TRUE)

existing <- extract_strata(sraster, existing)

oo <- sample_strat(sraster = sraster, nSamp = 20, existing = existing, include = TRUE, remove = TRUE)

test_that("Total outputs", {
  expect_equal(nrow(o), 20L)
  expect_equal(ncol(o), 4L)
  
  expect_equal(nrow(oo), 20L)
  expect_equal(ncol(oo), 4L)
  
})

test_that("messages", {
  
  expect_message(sample_strat(sraster = sraster, nSamp = 20, mindist = 200, access = access, buff_inner = 50, buff_outer = 200), "An access layer has been provided. An internal buffer of 50 m and an external buffer of 200 m have been applied.")
  expect_message(sample_strat(sraster = sraster, nSamp = 20, mindist = 200, access = access, buff_outer = 200), "An access layer has been provided. An external buffer of 200 m have been applied.")
  expect_message(sample_strat(sraster = sraster, nSamp = 20, existing = existing, include = TRUE, remove = TRUE), "'include = TRUE & remove = TRUE' - Stratum 1 overrepresented - 45 samples removed.")
  expect_message(sample_strat(sraster = sraster, nSamp = 200, existing = existing, include = TRUE, remove = TRUE), "Strata : 1 required no sample additions. Keeping all existing samples")
  expect_message(sample_strat(sraster = sraster, nSamp = 20, existing = data.frame(strata = c(1,2,3), x = c(1,2,3), y = c(1,2,3))),"'existing' column coordinate names are lowercase - converting to uppercase.")
  expect_message(sample_strat(sraster = sraster, nSamp = 25, access = access, buff_inner = 50, buff_outer = 200), "Buffered area contains 12454 available candidates. Sampling to reach 6 starting.")
  expect_message(sample_strat(sraster = sraster, nSamp = 20, filename = file.path(tempdir(), "temp.shp") , overwrite = TRUE), "Output samples written to disc.")
})

test_that("errors", {
  expect_error(sample_strat(sraster = sraster, nSamp = 20, include = TRUE),"'existing' must be provided when 'include = TRUE'.")
  expect_error(sample_strat(sraster = sraster, nSamp = 20, remove = TRUE),"'existing' must be provided when 'remove = TRUE'.")
  expect_error(sample_strat(sraster = sraster, nSamp = 20, existing = data.frame()),"'existing' must have an attribute named 'strata'. Consider using extract_strata().")
  expect_error(sample_strat(sraster = sraster, nSamp = 20, existing = data.frame(strata = c(1,2,3))),"'existing' must have columns named 'X' and 'Y'.")
  expect_error(sample_strat(sraster = sraster, nSamp = 15000, access = access, allocation = "equal", buff_inner = 50, buff_outer = 200), "Insufficient candidate samples within the buffered access extent. Consider altering buffer widths.")
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

