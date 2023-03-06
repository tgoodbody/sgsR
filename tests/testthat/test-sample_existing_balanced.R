# Test for correct output type

test_that("sample_existing_balanced should return an sf object", {
  expect_error(sample_existing_balanced(existing = existing_samples, nSamp = 10, algorithm = 1),"'algorithm' must be type character.")
  expect_error(sample_existing_balanced(existing = existing_samples, nSamp = 10, algorithm = "unknown"),"Unknown algorithm specified. Please use one of 'lpm2_kdtree', 'lcube', 'lcubestratified'.")
  expect_error(sample_existing_balanced(existing = existing_samples, nSamp = 10, p = "1"),"'p' must be type numeric.")
  expect_error(sample_existing_balanced(existing = existing_samples, nSamp = 10, p = 100),"'p' must have a length of 200.")
})


test_that("sample_existing_balanced should return an sf object", {
  samples <- sample_existing_balanced(existing = existing_samples, nSamp = 10)
  
  expect_equal(nrow(samples), 10L)
  expect_s3_class(samples, "sf")
})

# Test for correct output dimensions
test_that("sample_existing_balanced should return a sample with correct dimensions", {
  samples <- sample_existing_balanced(existing_samples, nSamp = 50)
  
  expect_equal(nrow(samples), 50)
})

# Test for balanced sampling
test_that("sample_existing_balanced should produce a balanced sample", {
  expect_error(sample_existing_balanced(existing, nSamp = 100, algorithm = "lcubestratified"),"'existing' must have a variable named 'strata' to use the 'lcubestratified' algorithm.")
  samples <- sample_existing_balanced(existing_samples, nSamp = 100, algorithm = "lcubestratified")
  
  # Check proportion of samples in each stratum
  counts <- table(existing_samples$strata)
  props <- counts / sum(counts)
  sample_counts <- table(samples$strata)
  sample_props <- sample_counts / sum(sample_counts)
  expect_equal(props, sample_props, tolerance = 0.1)
})

# Test for balanced sampling
test_that("lcube", {
  expect_s3_class(sample_existing_balanced(existing_samples, nSamp = 100, algorithm = "lcube"),"sf")
})