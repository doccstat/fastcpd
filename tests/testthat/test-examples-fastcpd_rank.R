testthat::test_that(
  "examples/fastcpd_rank.R", {
    source("examples/fastcpd_rank.R")
    testthat::expect_equal(result_mean@cp_set, c(332, 333))
    testthat::expect_equal(result_rank@cp_set, 200)
  }
)
