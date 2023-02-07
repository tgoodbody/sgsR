
test_that("Total outputs", {
  expect_error(calculate_pop(mraster = "A"), "'mraster' must be type SpatRaster.")
  expect_error(calculate_pop(mraster = mraster, PCA = 2), "'PCA' must be type logical.")
  expect_error(calculate_pop(mraster = mraster, matQ = 2), "'matQ' must be type logical.")
  expect_error(calculate_pop(mraster = mraster, matCov = 2), "'matCov' must be type logical.")
  
  expect_error(calculate_pop(mraster = mraster, nQuant = "A"), "'nQuant' must be type numeric.")
  expect_error(calculate_pop(mraster = mraster, matCov = TRUE, matQ = FALSE), "Covariance matrix creation requires quantile matrix. Set 'matQ = TRUE'.")
  
})

test_that("messages", {
  skip_on_cran()
  expect_equal(length(calculate_pop(mraster = mraster, matCov = FALSE)), 2L)
  
})

test_that("Total outputs", {
  skip_on_cran()
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


test_that("small raster", {
  skip_on_cran()
  
  mr1 <- mr <- rast(nrow = 50, ncol = 50)
  values(mr) <- runif(50*50)
  
  expect_equal(nrow(calculate_pop(mraster = mr)$values), 2500L)

  
  
})

