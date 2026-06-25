testthat::test_that(
  "examples/fastcpd_quantile_2.R", {
    source("examples/fastcpd_quantile_2.R")
    # lm reports spurious breaks at the outlier positions (50-52, 300-302)
    # plus the true break near 200.
    testthat::expect_equal(result_lm@cp_set, c(49, 50, 51, 52, 201, 299, 300, 301, 302))
    # quantile regression finds only the true change point.
    testthat::expect_equal(result_qr@cp_set, 201)
  }
)
