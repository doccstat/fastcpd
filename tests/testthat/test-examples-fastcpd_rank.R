testthat::test_that(
  "examples/fastcpd_rank.R", {
    source("examples/fastcpd_rank.R")
    testthat::expect_equal(result_mean@cp_set, c(332, 333))
    testthat::expect_equal(result_rank@cp_set, 200)
    testthat::expect_identical(result_rank@family, "rank")
    testthat::expect_equal(result_rank@data[[1]], x)
    testthat::expect_equal(rank_profile$estimate, result_rank@cp_set)
    centered_ranks <- rank(x) - (length(x) + 1) / 2
    bounds <- c(0, result_rank@cp_set, length(x))
    expected_se <- vapply(
      seq_len(length(bounds) - 1L),
      function(i) {
        segment <- centered_ranks[(bounds[i] + 1L):bounds[i + 1L]]
        stats::sd(segment) / sqrt(length(segment))
      },
      numeric(1)
    )
    testthat::expect_equal(rank_wald$se, expected_se)
  }
)
