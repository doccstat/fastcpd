testthat::test_that(
  "examples/fastcpd_kcp.R", {
    source("examples/fastcpd_kcp.R")
    testthat::expect_equal(result@cp_set, 200)
    testthat::expect_length(identical_result@cp_set, 0)
  }
)
