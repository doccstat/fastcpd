testthat::test_that(
  "examples/variance_arma.R", {
    source("examples/variance_arma.R")
    testthat::expect_equal(result$sigma2_bic, 100.25937)
  }
)
