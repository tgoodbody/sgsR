o <- sample_strat(sraster = sraster, nSamp = 20)
o1 <- sample_strat(sraster = sraster, nSamp = 20, allocation = "equal", details = TRUE)
o2 <- sample_strat(sraster = sraster, nSamp = 100, allocation = "manual", weights = c(0.2, 0.2, 0.2, 0.4), details = TRUE)
o3 <- sample_strat(sraster = sraster, mraster = mraster$zq90, nSamp = 100, allocation = "optim", weights = c(0.2, 0.2, 0.2, 0.4), details = TRUE, force = TRUE)

existing <- extract_strata(sraster, existing)
existingna <- extract_strata(sraster, existingna)

toSample2 <- toSample <- allocate_prop(sraster = sraster, nSamp = 100)

oo <- sample_strat(sraster = sraster, nSamp = 20, existing = existing, include = TRUE, remove = TRUE)

test_that("errors", {
  expect_error(sample_strat(sraster = "sraster", nSamp = 20), "'sraster' must be type SpatRaster.")
  expect_error(sample_strat(sraster = sraster, nSamp = "20"), "'nSamp' must be type numeric.")
  expect_error(sample_strat(sraster = sraster, nSamp = 20, method = 1), "'method' must be type character.")
  expect_error(sample_strat(sraster = sraster, nSamp = 20, method = "c"), "'method' must be one of 'random' or 'Queinnec'.")
  expect_error(sample_strat(sraster = sraster, nSamp = 20, mindist = "200"), "'mindist' must be type numeric.")
  expect_error(sample_strat(sraster = sraster, nSamp = 200, existing = existing, include = "TRUE", remove = TRUE), "'include' must be type logical.")
  expect_error(sample_strat(sraster = sraster, nSamp = 200, existing = existing, include = TRUE, remove = "TRUE"), "'remove' must be type logical.")
  expect_error(sample_strat(sraster = sraster, nSamp = 200, existing = existing, include = TRUE, remove = TRUE, force = "TRUE"), "'force' must be type logical.")
  expect_error(sample_strat(sraster = sraster, nSamp = 25, wrow = "A"), "'wrow' must be type numeric.")
  expect_error(sample_strat(sraster = sraster, nSamp = 25, wcol = "A"), "'wcol' must be type numeric.")
  expect_error(sample_strat(sraster = sraster, nSamp = 25, wrow = 2), "'wrow' must be an odd number.")
  expect_error(sample_strat(sraster = sraster, nSamp = 25, wcol = 2), "'wcol' must be an odd number.")
  expect_error(sample_strat(sraster = sraster, nSamp = 25, plot = "TRUE"), "'plot' must be type logical.")
  expect_error(sample_strat(sraster = sraster, nSamp = 25, details = "TRUE"), "'details' must be type logical.")
  expect_error(sample_strat(sraster = sraster, nSamp = 25, existing = "A"), "'existing' must be a data.frame or sf object.")
  expect_error(sample_strat(sraster = sraster, nSamp = 25, existing = access), "'existing' geometry type must be 'sfc_POINT'.")

  expect_error(sample_strat(sraster = sraster, nSamp = 20, include = TRUE), "'existing' must be provided when 'include = TRUE'.")
  expect_error(sample_strat(sraster = sraster, nSamp = 20, remove = TRUE), "'existing' must be provided when 'remove = TRUE'.")
  expect_error(sample_strat(sraster = sraster, nSamp = 20, existing = data.frame()), "'existing' must have an attribute named 'strata'. Consider using extract_strata().")
  expect_error(sample_strat(sraster = sraster, nSamp = 20, existing = data.frame(strata = c(1, 2, 3))), "'existing' must have columns named 'X' and 'Y'.")
  expect_error(sample_strat(sraster = sraster, nSamp = 15000, access = access, allocation = "equal", buff_inner = 50, buff_outer = 200), "Insufficient candidate samples within the buffered access extent. Consider altering buffer widths.")
  expect_error(sample_strat(sraster = sraster, nSamp = 20, filename = file.path(tempdir(), "temp.shp"), overwrite = "TRUE"), "'overwrite' must be type logical.")
})


test_that("Total outputs", {
  expect_equal(nrow(o), 20L)
  expect_equal(ncol(o), 4L)

  expect_equal(nrow(oo), 20L)
  expect_equal(ncol(oo), 4L)

  #--- categorical ---#
  expect_message(sample_strat(sraster = x, nSamp = 200), "'sraster' has factor values. Converting to allow mapping.")
  expect_equal(ncol(sample_strat(sraster = x, nSamp = 200, access = access, buff_outer = 100, plot = TRUE)), 5L)
  expect_equal(ncol(sample_strat(sraster = x, nSamp = 200, plot = TRUE)), 5L)
})

test_that("messages", {
  expect_message(sample_strat(sraster = sraster, nSamp = 20), "Using 'Queinnec' sampling method.")
  expect_message(sample_strat(sraster = sraster, nSamp = 20, method = "random"), "Using 'random' sampling method. Ignoring 'existing', 'include', 'remove' if provided.")
  expect_message(sample_strat(sraster = sraster, nSamp = 20, mindist = 200, access = access, buff_inner = 50, buff_outer = 200), "An access layer has been provided. An internal buffer of 50 m and an external buffer of 200 m have been applied.")
  expect_message(sample_strat(sraster = sraster, nSamp = 20, mindist = 200, access = access, buff_outer = 200), "An access layer has been provided. An external buffer of 200 m have been applied.")
  expect_message(sample_strat(sraster = sraster, nSamp = 20, existing = existing, include = TRUE, remove = TRUE), "'include = TRUE & remove = TRUE' - Stratum 1 overrepresented - 45 samples removed.")
  expect_message(sample_strat(sraster = sraster, nSamp = 20, existing = existing, include = TRUE, remove = FALSE), "'include = TRUE & remove = FALSE' - Stratum 1 overrepresented by 45 samples but have not been removed. Expect a higher total 'nSamp' in output.")
  expect_message(sample_strat(sraster = sraster, nSamp = 200, existing = existing, include = TRUE, remove = TRUE), "Strata : 1 required no sample additions. Keeping all existing samples")
  expect_message(sample_strat(sraster = sraster, nSamp = 20, existing = data.frame(strata = c(1, 2, 3), x = c(1, 2, 3), y = c(1, 2, 3))), "'existing' column coordinate names are lowercase - converting to uppercase.")
  expect_message(sample_strat(sraster = sraster, nSamp = 20, existing = existingna), "Implementing proportional allocation of samples.")

  expect_message(sample_strat(sraster = sraster, nSamp = 25, access = access, buff_inner = 50, buff_outer = 200), "Buffered area contains 12454 available candidates. Sampling to reach 6 starting.")
  expect_message(sample_strat(sraster = sraster, nSamp = 20, filename = file.path(tempdir(), "temp.shp"), overwrite = TRUE), "Output samples written to disc.")
})

test_that("Test equal", {
  expect_message(sample_strat(sraster = sraster, mraster = mraster$zq90, nSamp = 10, allocation = "equal", weights = c(0.2, 0.2, 0.2, 0.4)), "'weights' was specified but 'allocation = equal' - did you mean to use 'allocation = manual'?")
  expect_message(sample_strat(sraster = sraster, mraster = mraster$zq90, nSamp = 10, allocation = "equal", weights = c(0.2, 0.2, 0.2, 0.4)), "'mraster' was specified but 'allocation = equal' - did you mean to use 'allocation = optim'?")
  expect_equal(nrow(o1$sampleDist), 4L)
  expect_equal(sum(o1$sampleDist$total), 80L)
  expect_equal(unique(o1$samples$type), "new")

  expect_equal(sum(o2$sampleDist$total), 100L)

  expect_type(o2, "list")
  expect_s3_class(o1$samples, "sf")

  expect_equal(nrow(sample_strat(sraster = sraster, allocation = "equal", nSamp = 5, method = "random")), 20L)
})

test_that("Test manual", {
  expect_message(sample_strat(sraster = sraster, mraster = mraster$zq90, nSamp = 10, allocation = "manual", weights = c(0.2, 0.2, 0.2, 0.4)), "'mraster' was specified but 'allocation = manual' - did you mean to use 'allocation = optim'?")
  expect_error(sample_strat(sraster = sraster, nSamp = 10, allocation = "manual"), "'weights' must be defined if 'allocation = manual'.")

  expect_equal(o2$sampleDist$total[1], 20L)

  expect_type(o2, "list")
  expect_s3_class(o2$samples, "sf")

  expect_equal(nrow(sample_strat(sraster = sraster, allocation = "manual", weights = c(0.2, 0.2, 0.2, 0.4), nSamp = 10, method = "random")), 10L)
})

test_that("Test optim", {
  expect_message(sample_strat(sraster = sraster, mraster = mraster$zq90, nSamp = 10, allocation = "optim", weights = c(0.2, 0.2, 0.2, 0.4), details = TRUE), "'weights' was specified but 'allocation = optim' - did you mean to use 'allocation = manual'?")
  expect_error(sample_strat(sraster = sraster, nSamp = 20, allocation = "optim", mraster = "m"), "'mraster' must be type SpatRaster.")
  expect_message(sample_strat(sraster = sraster, mraster = mraster$zq90, nSamp = 100, allocation = "optim", weights = c(0.2, 0.2, 0.2, 0.4), details = TRUE), "nSamp of 100 is not perfectly divisible based on strata distribution. nSamp of 99 will be returned. Use 'force = TRUE' to brute force to 100.")
  expect_equal(sum(o3$sampleDist$total), 100L)

  expect_equal(sum(o2$sampleDist$total), 100L)

  expect_type(o3, "list")
  expect_s3_class(o3$samples, "sf")

  expect_equal(nrow(sample_strat(sraster = sraster, allocation = "optim", mraster = mraster$zq90, nSamp = 20, method = "random")), 20L)
})

test_that("existing", {
  expect_message(allocate_existing(toSample = toSample, existing = existingna), "16 samples in 'existing' are located where strata values are NA. Expect 16 additional samples in output.")

  toSample2$strata <- c(5, 6, 7, 8)

  expect_error(allocate_existing(toSample = toSample2, existing = existingna), "'existing' does not contain matching strata to those in 'sraster'. Check strata in both data sets & consider using extract_strata().")
})

test_that("force", {
  expect_equal(sum(allocate_force(toSample = toSample, nSamp = 100, diff = 1)$total), 99L)
  expect_equal(sum(allocate_force(toSample = toSample, nSamp = 100, diff = -1)$total), 101L)
})


test_that("category column", {
  xx <- strat_map(c(x, sraster))

  expect_equal(ncol(sample_strat(xx, nSamp = 1000, method = "random")), 3)
})
