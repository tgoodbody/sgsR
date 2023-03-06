set.seed(2022)
o <- sample_balanced(mraster = mraster, nSamp = 50)

xrna <- rast(nrow = 1, ncol = 1, crs = NA)
values(xrna) <- 1

test_that("inputs", {
  expect_error(sample_balanced(mraster = "A", nSamp = 4), "'mraster' must be type SpatRaster.")
  expect_error(sample_balanced(mraster = mraster, nSamp = "4"), "'nSamp' must be type numeric.")
  expect_error(sample_balanced(mraster = mraster, nSamp = 4, plot = "TRUE"), "'plot' must be type logical.")
  expect_error(sample_balanced(mraster = mraster, nSamp = 4, algorithm = 3), "'algorithm' must be type character.")
  expect_error(sample_balanced(mraster = mraster, nSamp = 4, algorithm = "TRUE"), "Unknown algorithm specified. Please use one of 'lpm2_kdtree', 'lcube', 'lcubestratified'.")
  expect_error(sample_balanced(mraster = mraster, nSamp = 4, p = "p"), "'p' must be type numeric.")
  expect_error(sample_balanced(mraster = mraster, nSamp = 4, p = c(1, 2, 3)), "'p' have a length of 91195.")
})



test_that("Total outputs", {
  expect_equal(nrow(o), 50L)
  expect_equal(ncol(o), 1L)
  expect_s3_class(o, "sf")
  
  sample_clhs(mraster = mraster, nSamp = 320, access = access, buff_inner = 50, buff_outer = 200, plot = TRUE)
})

test_that("errors", {
  expect_error(sample_clhs(mraster = mraster, nSamp = 320, cost = "A"), "No layer named 'A' exists in 'mraster'.")
  expect_error(sample_balanced(mraster = mraster, nSamp = 50, algorithm = "lcubestratified"), "'mraster' must have a variable named 'strata' to use the 'lcubestratified' algorithm")
})

test_that("Messages", {
  skip_on_cran()
  expect_message(sample_balanced(mraster = mraster, nSamp = 50, access = access, buff_inner = 50, buff_outer = 200, plot = TRUE), "An access layer has been provided. An internal buffer of 50 m and an external buffer of 200 m have been applied.")
  expect_message(sample_balanced(mraster = mraster, nSamp = 50, access = access, buff_outer = 200), "An access layer has been provided. An external buffer of 200 m have been applied.")
})


test_that("lcube", {
  skip_on_cran()
  expect_equal(nrow(sample_balanced(mraster = mraster, nSamp = 50, algorithm = "lcube")), 50L)
})

test_that("lcubestratified", {
  skip_on_cran()

  set.seed(2022)
  mrs <- c(mraster, sraster)

  o1 <- sample_balanced(mraster = mrs, nSamp = 50, algorithm = "lcubestratified")

  expect_equal(nrow(o1), 50L)
})
