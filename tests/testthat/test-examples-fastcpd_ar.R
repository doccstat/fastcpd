testthat::test_that(
  "examples/fastcpd_ar.R", {
    source("examples/fastcpd_ar.R")
    testthat::expect_equal(result@cp_set, 614)
  }
)
