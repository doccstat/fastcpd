testthat::test_that(
  "examples/fastcpd_poisson.R", {
    testthat::skip_if_not_installed("mvtnorm")
    source("examples/fastcpd_poisson.R")
    testthat::expect_equal(result@cp_set, c(312, 716, 1020, 1102, 1301))
  }
)
