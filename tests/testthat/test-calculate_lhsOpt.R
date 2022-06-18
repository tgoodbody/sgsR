test_that("Total outputs", {
  skip_on_cran()
  
  set.seed(2022)
  #--- supply quantile and covariance matrices ---#
  mats <- calculate_pop(mraster = mraster, PCA = TRUE)
  
  #--- supply quantile and covariance matrices ---#
  expect_message(m <- calculate_lhsOpt(mats = mats),"Your optimum estimated sample size based on KL divergence is: 40")
  
  expect_equal(nrow(m), 10L)
  expect_equal(ncol(m), 7L)
  expect_equal(names(m),c("n", "mean_dist", "sd_dist", "min_S", "max_S", "mean_KL", "sd_KL"))

})

