test_that("Total outputs", {
  set.seed(2022)
  #--- supply quantile and covariance matrices ---#
  mats <- calculate_pop(mraster = mraster, PCA = TRUE)
  
  expect_equal(nrow(mats$values), 91195L)
  expect_equal(nrow(mats$matQ), 11L)
  expect_equal(nrow(mats$matCov), 10L)
  expect_equal(ncol(mats$values), 3L)
  expect_equal(length(mats), 4L)
  expect_equal(names(mats),c("values","pcaLoad","matQ","matCov"))
  
  expect_named(mats$values,c("zq90","pzabove2","zsd"))
  expect_equal(sum(mats$matQ[,1]),169.399996)
  expect_equal(sum(mats$matQ[,2]),550)
})
