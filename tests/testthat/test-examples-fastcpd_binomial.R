testthat::test_that(
  "examples/fastcpd_binomial.R", {
    source("examples/fastcpd_binomial.R")
    testthat::expect_equal(result@cp_set, 126)
  }
)
