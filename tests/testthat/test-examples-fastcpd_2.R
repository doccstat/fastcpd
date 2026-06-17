testthat::test_that(
  "examples/fastcpd_2.R", {
    testthat::skip_if(Sys.getenv("R_COVR") == "true", "slow custom Huber example: skip during coverage")
    testthat::skip_if_not_installed("mvtnorm")
    testthat::skip_if_not_installed("stats")

    source("examples/fastcpd_2.R")
    testthat::expect_equal(huber_regression_result@cp_set, c(418, 726))
  }
)
