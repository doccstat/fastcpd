testthat::test_that("plot accepts change-point-only fits", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::expect_no_error(source("examples/plot_cp_only.R"))
  testthat::expect_equal(plot_cp_only_fit@cp_set, 50)
  testthat::expect_equal(plot_cp_only_fit@cp_set, plot_detail_fit@cp_set)
  testthat::expect_length(plot_cp_only_fit@cost_values, 0)
  testthat::expect_equal(length(plot_cp_only_fit@residuals), 0)
  testthat::expect_equal(nrow(plot_detail_fit@residuals), length(plot_data))
})
