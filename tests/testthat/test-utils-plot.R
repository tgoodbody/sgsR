o <- plot_scatter(mraster = mraster, existing = existing, samp = 0.0001)
o1 <- plot_scatter(mraster = mraster, existing = existing, reverse = TRUE, samp = 0.0001)

test_that("scatter errors", {
  expect_error(plot_scatter(mraster = "mraster", existing = "a"), "'mraster' must be type SpatRaster.")
  expect_error(plot_scatter(mraster = mraster, existing = "a"), "'existing' must be an sf object.")
  expect_error(plot_scatter(mraster = mraster, existing = existing, reverse = 1), "'reverse' must be type logical.")
  expect_error(plot_scatter(mraster = mraster, existing = existing, samp = "a"), "'samp' must be type numeric.")

  expect_error(plot_scatter(mraster = mraster[[1]], existing = existing), "Only 1 layer in `mraster` when 2 are needed.")
  expect_error(plot_scatter(mraster = mraster[[1:2]], existing = existing, reverse = "1"), "'reverse' must be type logical.")
  expect_error(plot_scatter(mraster = mraster, existing = existing, samp = 1.1), "'samp' must be > 0 <= 1.")
})

test_that("scatter messages", {
  expect_message(plot_scatter(mraster = mraster, existing = existing), "More than 2 layers in `mraster`. Only first 2 layers will be used.")
  expect_s3_class(o, "gg")
  
  get_labs <- function(x) x$labels
  if ("get_labs" %in% getNamespaceExports("ggplot2")) {
    get_labs <- ggplot2::get_labs
  }
  
  expect_equal(get_labs(o)$x, "zq90")
  expect_equal(get_labs(o1)$x, "pzabove2")
})
