testthat::test_that(
  "ggplot2 is not installed selection 2", {
    # TODO(doccstat): I can not find a way to cover the `if` branch in
    # `plot.fastcpd` when `ggplot2` is not installed. `stub` from `mockery`
    # does not work.
    testthat::skip("This test is intended to be run manually.")
    match_call <- ""
    class(match_call) <- "language"
    thetas <- data.frame(matrix(c(1, 2), 1, 2))
    names(thetas) <- c("segment 1", "segment 2")

    class_instance <- methods::new(
      Class = "fastcpd",
      call = match_call,
      data = data.frame(y = c(1, 1, 2, 2), x = rep(1, 4)),
      family = "lasso",
      cp_set = 3,
      cost_values = c(9.5, 10.5),
      residuals = rep(0, 4),
      thetas = thetas,
      cp_only = FALSE
    )

    mockery::stub(plot, "utils::menu", function(...) 2)
    mockery::stub(plot, "base::requireNamespace", function(...) FALSE)

    testthat::expect_message(
      plot(class_instance),
      "ggplot2 is not installed. No plot is made.\n"
    )
  }
)

testthat::test_that(
  "ggplot2 is not installed selection 1", {
    match_call <- ""
    class(match_call) <- "language"
    thetas <- data.frame(matrix(c(1, 2), 1, 2))
    names(thetas) <- c("segment 1", "segment 2")

    class_instance <- methods::new(
      Class = "fastcpd",
      call = match_call,
      data = data.frame(y = c(1, 1, 2, 2), x = rep(1, 4)),
      family = "lasso",
      cp_set = 3,
      cost_values = c(9.5, 10.5),
      residuals = rep(0, 4),
      thetas = thetas,
      cp_only = FALSE
    )

    # TODO(doccstat): I can not find a way to cover the `if` branch in
    # `plot.fastcpd` when `ggplot2` is not installed. `menu` with selection 1
    # should be mocked here.
    testthat::expect_no_error(plot(class_instance))

    testthat::expect_equal(
      capture.output(print(class_instance)),
      c(
        "",
        "Change points:",
        "[1] 3"
      )
    )

    testthat::expect_equal(
      capture.output(show(class_instance)),
      c(
        "",
        "A fastcpd object.",
        "Available methods to evaluate the object are:",
        "plot, print, show, summary",
        "",
        "",
        "Change points:",
        "[1] 3"
      )
    )

    testthat::expect_equal(
      capture.output(summary(class_instance)),
      c(
        "",
        "Call:",
        "structure(\"\", class = \"language\")",
        "",
        "Change points:",
        "3 ",
        "",
        "Cost values:",
        "9.5 10.5 ",
        "",
        "Parameters:",
        "1 x 2 sparse Matrix of class \"dgCMatrix\"",
        "     segment 1 segment 2",
        "[1,]         1         2"
      )
    )

    class_instance@family <- "custom"

    testthat::expect_equal(
      capture.output(summary(class_instance)),
      c(
        "",
        "Call:",
        "structure(\"\", class = \"language\")",
        "",
        "Change points:",
        "3 ",
        "",
        "Parameters:",
        "  segment 1 segment 2",
        "1         1         2"
      )
    )
  }
)

testthat::test_that("output methods without change points", {
  match_call <- ""
  class(match_call) <- "language"

  class_instance <- methods::new(
    Class = "fastcpd",
    call = match_call,
    data = data.frame(matrix(NA, 0, 0)),
    family = "lasso",
    cp_set = numeric(0),
    cost_values = numeric(0),
    residuals = numeric(0),
    thetas = data.frame(matrix(NA, 0, 0)),
    cp_only = TRUE
  )

  testthat::expect_equal(
    capture.output(print(class_instance)),
    c(
      "",
      "No change points found"
    )
  )

  testthat::expect_equal(
    capture.output(summary(class_instance)),
    c(
      "",
      "Call:",
      "structure(\"\", class = \"language\")",
      "",
      "No change points found"
    )
  )
})
