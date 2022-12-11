xna <- rast(ncol = 4, nrow = 4)
values(xna) <- NA


strata_m <- terra::mask(sraster,
                        mask = sraster,
                        maskvalues = 1,
                        inverse = TRUE)

add_strata <- strat_rule1(n = 200,s = 1,i = 1, strat_mask = strata_m, add_strata = data.frame(), extraCols = NULL, mindist = 100)


test_that("strat rule 1", {
  
  expect_equal(add_strata$add_strata$strata[1], 1L)
  expect_equal(strat_rule1(n = 200,s = 1,i = 2, strat_mask = strata_m, add_strata = add_strata$add_strata, extraCols = NULL, mindist = 100)$add_strata$strata[1], 1L)
  expect_equal(nrow(strat_rule1(n = 200,s = 1,i = 2, strat_mask = xna, add_strata = data.frame(), extraCols = NULL))$add_strata, NULL)
})


test_that("strat rule 2", {

  add_strata_NA <- strat_rule2(n = 200, s = 1, nCount = 0, strata_m = strata_m, add_strata = data.frame(), extraCols = NULL, mindist = NULL)
  add_strata_nonNA <- strat_rule2(n = 200, s = 1, nCount = add_strata$nCount, strata_m = strata_m, add_strata = add_strata$add_strata, extraCols = NULL, mindist = NULL)
  
  expect_equal(nrow(add_strata_NA), 200L)
  
  expect_equal(nrow(strat_rule2(n = 200, s = 1, nCount = 0, strata_m = xna, add_strata = data.frame(), extraCols = NULL, mindist = NULL)), 0L)

})