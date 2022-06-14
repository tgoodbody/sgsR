test_that("Total outputs", {
  expect_equal(nrow(mat$values), 91195L)
  expect_equal(nrow(mat$matQ), 11L)
  expect_equal(nrow(mat$matCov), 10L)
  expect_equal(ncol(mat$values), 3L)
  
  expect_named(mat$values,c("zq90","pzabove2","zsd"))
  expect_equal(sum(mat$matQ[,1]),174.899996)
  expect_equal(sum(mat$matQ[,2]),548.35003)
})
