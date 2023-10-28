testthat::test_that(
  "man/examples/fastcpd_lm.R", {
    source("man/examples/fastcpd_lm.R")
    testthat::expect_equal(result_internal_variance@cp_set, c(100, 201))
  }
)
