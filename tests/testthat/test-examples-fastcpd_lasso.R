testthat::test_that(
  "examples/fastcpd_lasso.R", {
    testthat::skip_if_not_installed("mvtnorm")
    source("examples/fastcpd_lasso.R")
    testthat::expect_equal(result@cp_set, c(80, 200, 320))
  }
)
