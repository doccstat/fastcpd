testthat::test_that("plot uses the first aligned residual response", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::expect_no_error(source("examples/plot_residuals.R"))
  testthat::expect_equal(nrow(plot_residual_fit@residuals), length(plot_x))
  testthat::expect_equal(ncol(plot_multiple_residual_fit@residuals), 2)
  testthat::expect_equal(plot_regular_residuals, plot_residual_fit@residuals[, 1])
  testthat::expect_identical(plot_first_residuals, plot_regular_residuals)
  for (residuals in list(
    plot_missing_residuals, plot_short_residuals, plot_zero_column_residuals
  )) {
    testthat::expect_length(residuals, length(plot_x))
    testthat::expect_true(all(is.na(residuals)))
  }
})
