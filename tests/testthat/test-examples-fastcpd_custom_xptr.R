testthat::test_that(
  "examples/fastcpd_custom_xptr.txt", {
    testthat::skip_on_cran()
    testthat::skip_if_not_installed("Rcpp")
    testthat::skip_if_not_installed("RcppArmadillo")

    examples_xptr <- readLines("examples/fastcpd_custom_xptr.txt")
    tryCatch(
      source(textConnection(paste(
        examples_xptr[seq_len(length(examples_xptr) - 2) + 1],
        collapse = "\n"
      ))),
      error = function(e) {
        testthat::skip(paste("C++ compilation failed:", conditionMessage(e)))
      }
    )

    testthat::expect_equal(result_xptr@cp_set, 500)
    testthat::expect_equal(
      unname(unlist(result_xptr@cost_values)),
      c(255.487731877, 278.914161918)
    )
    testthat::expect_identical(result_r@cp_set, result_builtin@cp_set)
    testthat::expect_identical(result_xptr@cp_set, result_builtin@cp_set)
  }
)
