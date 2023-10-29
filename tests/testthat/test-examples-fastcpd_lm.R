testthat::test_that(
  "examples/fastcpd_lm.R", {
    source("examples/fastcpd_lm.R")
    testthat::expect_equal(result@cp_set, c(97, 201))
  }
)
