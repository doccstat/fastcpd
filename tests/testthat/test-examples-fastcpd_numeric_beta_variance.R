testthat::test_that(
  "numeric Gaussian penalties skip unused automatic variance estimation", {
    source("examples/fastcpd_numeric_beta_variance.R")

    testthat::expect_s4_class(numeric_beta_lm_result, "fastcpd")
    testthat::expect_equal(
      numeric_beta_lm_result@cp_set,
      numeric_beta_lm_explicit_variance@cp_set
    )
  }
)
