testthat::test_that(
  "examples/fastcpd_var.R", {
    testthat::skip_on_cran()
    source("examples/fastcpd_var.R")
    testthat::expect_equal(result@cp_set, 204)
  }
)
