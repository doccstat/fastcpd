testthat::test_that(
  "examples/fastcpd_quantile.txt", {
    source("examples/fastcpd_quantile.txt")
    testthat::expect_equal(result@cp_set, 100)
    testthat::expect_equal(result@thetas, data.frame(
      "segment 1" = c(0.9882078, -1.0449301),
      "segment 2" = c(-1.1610592, 0.9042203),
      check.names = FALSE
    ), tolerance = 1e-5)
  }
)
