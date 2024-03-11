testthat::test_that(
  "examples/variance_median.R", {
    source("examples/variance_median.R")
    testthat::expect_equal(sigma2, 5473705.3)
  }
)
