testthat::test_that(
  "examples/fastcpd_poisson.txt", {
    testthat::skip_if_not_installed("ggplot2")
    testthat::skip_if_not_installed("mvtnorm")

    examples_poisson <- readLines("examples/fastcpd_poisson.txt")
    source(textConnection(paste(
      examples_poisson[seq_len(length(examples_poisson) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, c(506, 838, 1003))
    testthat::expect_equal(result@thetas, data.frame(
      "segment 1" = c(1.015468142, 0.276378274, -1.049326186),
      "segment 2" = c(0.656870530, -0.213138662, -0.594279490),
      "segment 3" = c(1.037186063, 0.264881329, -0.980155436),
      "segment 4" = c(1.445192824, 0.991007927, -1.435463755),
      check.names = FALSE
    ))
  }
)
