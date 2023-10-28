testthat::test_that(
  "man/examples/fastcpd_binomial.R", {
    source("man/examples/fastcpd_binomial.R")
    testthat::expect_equal(result_ts@cp_set, 126)
  }
)
