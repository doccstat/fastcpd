testthat::test_that(
  "examples/fastcpd_lm_1.R", {
    testthat::skip_if_not_installed("ggplot2")
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/fastcpd_lm_1.R")
    testthat::expect_equal(result_lm@cp_set, c(97, 201))
    testthat::expect_equal(result_lm@thetas, data.frame(
      "segment 1" = c(0.74291290, 3.69465275, -1.24746871, 0.09579985),
      "segment 2" = c(-0.615304925, -0.503494764, 2.252235245, -1.987512620),
      "segment 3" = c(0.873347311, 0.322286842, 1.018845519, 2.276134007),
      check.names = FALSE
    ))
  }
)
