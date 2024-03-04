testthat::test_that(
  "examples/variance_arma.R", {
    source("examples/variance_arma.R")
    testthat::expect_equal(sigma2, 102.892884)
  }
)
