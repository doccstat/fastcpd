# These tests target `FitGlm` (`src/fastcpd_glm.h`), the from-scratch IRLS
# GLM fitter that replaced the vendored `ref_fastglm_*`/`ref_r_family.c` code.
# IRLS for a canonical-link GLM converges to the unique MLE, so on data with
# no real change point and a penalty large enough to forbid splitting,
# `fastcpd`'s single-segment estimate should match `glm()`/`lm()` to high
# precision -- this directly validates `FitGlm` against R's reference
# implementation rather than only checking `fastcpd`'s end-to-end output.

testthat::test_that(
  "FitGlm matches glm() MLE for binomial family on a single segment", {
    set.seed(1)
    n <- 200
    p <- 2
    x <- matrix(rnorm(n * p), ncol = p)
    theta_true <- c(0.6, -0.9)
    y <- rbinom(n, 1, 1 / (1 + exp(-x %*% theta_true)))

    result <- suppressWarnings(
      fastcpd.binomial(cbind(y, x), beta = 1e6)
    )
    glm_fit <- stats::glm(
      y ~ . - 1, data = data.frame(y = y, x = x), family = stats::binomial()
    )

    testthat::expect_length(result@cp_set, 0)
    testthat::expect_equal(
      unname(unlist(result@thetas)),
      unname(stats::coef(glm_fit)),
      tolerance = 1e-4
    )
  }
)

testthat::test_that(
  "FitGlm matches glm() MLE for poisson family on a single segment", {
    set.seed(1)
    n <- 200
    p <- 2
    x <- matrix(rnorm(n * p, 0, 0.3), ncol = p)
    theta_true <- c(0.4, -0.3)
    y <- rpois(n, exp(x %*% theta_true))

    result <- suppressWarnings(
      fastcpd.poisson(cbind(y, x), beta = 1e6)
    )
    glm_fit <- stats::glm(
      y ~ . - 1, data = data.frame(y = y, x = x), family = stats::poisson()
    )

    testthat::expect_length(result@cp_set, 0)
    testthat::expect_equal(
      unname(unlist(result@thetas)),
      unname(stats::coef(glm_fit)),
      tolerance = 1e-4
    )
  }
)

testthat::test_that(
  "FitGlm matches lm() MLE for gaussian family on a single segment", {
    set.seed(1)
    n <- 200
    p <- 2
    x <- matrix(rnorm(n * p), ncol = p)
    theta_true <- c(1.2, -0.7)
    y <- x %*% theta_true + rnorm(n, 0, 0.5)

    result <- suppressWarnings(
      fastcpd.lm(cbind(y, x), beta = 1e6)
    )
    lm_fit <- stats::lm(y ~ . - 1, data = data.frame(y = y, x = x))

    testthat::expect_length(result@cp_set, 0)
    testthat::expect_equal(
      unname(unlist(result@thetas)),
      unname(stats::coef(lm_fit)),
      tolerance = 1e-4
    )
  }
)
