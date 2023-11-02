testthat::test_that(
  "ggplot2 is not installed mock", {
    testthat::skip_if_not_installed("mockthat")
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

    testthat::expect_error(
      mockthat::with_mock(
        `require_namespace` = function(...) FALSE,
        `utils_menu` = function(...) 2,
        plot(class_instance)
      ),
      "ggplot2 is not installed. No plot is made."
    )

    testthat::expect_no_error(
      mockthat::with_mock(
        `require_namespace` = function(...) FALSE,
        `utils_menu` = function(...) 1,
        `install_packages` = function(...) TRUE,
        plot(class_instance)
      )
    )

    testthat::expect_error(
      mockthat::with_mock(
        `require_namespace` = function(...) FALSE,
        `utils_menu` = function(...) 1,
        `install_packages` = function(...) {
          stop("ggplot2 could not be installed.")
        },
        plot(class_instance)
      ),
      "ggplot2 could not be installed."
    )
  }
)

testthat::test_that(
  "utility functions output test", {
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

testthat::test_that("mean change p > 1", {
  match_call <- ""
  class(match_call) <- "language"

  class_instance <- methods::new(
    Class = "fastcpd",
    call = match_call,
    data = data.frame(matrix(NA, 0, 2)),
    family = "mean",
    cp_set = numeric(0),
    cost_values = numeric(0),
    residuals = numeric(0),
    thetas = data.frame(matrix(NA, 0, 0)),
    cp_only = FALSE
  )

  testthat::expect_error(
    plot(class_instance),
    "Can not plot mean change points with p > 1."
  )
})
