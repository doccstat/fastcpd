testthat::test_that(
  "examples/variance_median.R", {
    testthat::skip_on_cran()
    source("examples/variance_median.R")
    testthat::expect_equal(sigma2, 5803645.3)
  }
)
