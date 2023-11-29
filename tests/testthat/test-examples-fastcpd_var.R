testthat::test_that(
  "examples/fastcpd_var.R", {
    source("examples/fastcpd_var.R")
    testthat::expect_equal(result@cp_set, 519)
  }
)
