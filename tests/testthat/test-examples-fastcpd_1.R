testthat::test_that(
  "examples/fastcpd_1.R", {
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/fastcpd_1.R")
    testthat::expect_equal(result_mlm@cp_set, 125)
  }
)
