# These tests validate the compiled (`Rcpp::XPtr`) custom cost function path
# added alongside the long-standing R-closure `family = "custom"` mechanism.
# `fastcpd_xptr.h` detects `EXTPTRSXP` SEXPs passed through the existing
# `cost` / `cost_gradient` / `cost_hessian` arguments (vs. R closures) and
# unwraps tagged function pointers directly -- bypassing R-call overhead.
# Compiled and R-closure paths compute the *same* mathematical cost, so a
# faithful implementation must produce numerically identical PELT results.

testthat::skip_on_cran()
testthat::skip_if_not_installed("Rcpp")

xptr_helpers_env <- new.env()
xptr_helpers <- tryCatch(
  {
    Rcpp::sourceCpp(
      code = '
      // [[Rcpp::depends(RcppArmadillo)]]
      #include <RcppArmadillo.h>

      double mean_nll_pelt(arma::mat const& data) {
        arma::rowvec const seg_mean = arma::mean(data, 0);
        arma::mat const centered = data.each_row() - seg_mean;
        return arma::accu(arma::square(centered)) / 2.0;
      }
      typedef double (*CostPeltFnPtr)(arma::mat const&);

      // [[Rcpp::export]]
      SEXP fastcpd_test_make_cost_pelt_xptr() {
        Rcpp::XPtr<CostPeltFnPtr> xptr(
            new CostPeltFnPtr(&mean_nll_pelt), true,
            Rcpp::wrap("fastcpd_cost_pelt"));
        return xptr;
      }

      double mean_nll_sen(arma::mat const& data, arma::colvec const& theta) {
        arma::vec const resid = data.col(0) - theta(0);
        return arma::accu(arma::square(resid)) / 2.0;
      }
      typedef double (*CostSenFnPtr)(arma::mat const&, arma::colvec const&);

      // [[Rcpp::export]]
      SEXP fastcpd_test_make_cost_sen_xptr() {
        Rcpp::XPtr<CostSenFnPtr> xptr(
            new CostSenFnPtr(&mean_nll_sen), true,
            Rcpp::wrap("fastcpd_cost_sen"));
        return xptr;
      }

      arma::colvec mean_gradient(arma::mat const& data, arma::colvec const& theta) {
        arma::vec const resid = data.col(0) - theta(0);
        arma::colvec gradient(1);
        gradient(0) = -arma::accu(resid);
        return gradient;
      }
      typedef arma::colvec (*CostGradientFnPtr)(arma::mat const&, arma::colvec const&);

      // [[Rcpp::export]]
      SEXP fastcpd_test_make_cost_gradient_xptr() {
        Rcpp::XPtr<CostGradientFnPtr> xptr(
            new CostGradientFnPtr(&mean_gradient), true,
            Rcpp::wrap("fastcpd_cost_gradient"));
        return xptr;
      }

      arma::mat mean_hessian(arma::mat const& data, arma::colvec const& theta) {
        arma::mat hessian(1, 1);
        hessian(0, 0) = static_cast<double>(data.n_rows);
        return hessian;
      }
      typedef arma::mat (*CostHessianFnPtr)(arma::mat const&, arma::colvec const&);

      // [[Rcpp::export]]
      SEXP fastcpd_test_make_cost_hessian_xptr() {
        Rcpp::XPtr<CostHessianFnPtr> xptr(
            new CostHessianFnPtr(&mean_hessian), true,
            Rcpp::wrap("fastcpd_cost_hessian"));
        return xptr;
      }
    ',
      env = xptr_helpers_env
    )
    xptr_helpers_env
  },
  error = function(e) NULL
)
testthat::skip_if(is.null(xptr_helpers), "no C++ compiler available for test fixtures")

fastcpd_test_make_cost_pelt_xptr <- xptr_helpers_env$fastcpd_test_make_cost_pelt_xptr
fastcpd_test_make_cost_sen_xptr <- xptr_helpers_env$fastcpd_test_make_cost_sen_xptr
fastcpd_test_make_cost_gradient_xptr <- xptr_helpers_env$fastcpd_test_make_cost_gradient_xptr
fastcpd_test_make_cost_hessian_xptr <- xptr_helpers_env$fastcpd_test_make_cost_hessian_xptr

set.seed(1)
data <- matrix(c(rnorm(150, mean = 0), rnorm(150, mean = 5)))

testthat::test_that(
  "PELT-style (one-argument) compiled cost matches its R-closure equivalent", {
    cost_pelt_r <- function(data) {
      sum((data - mean(data))^2) / 2
    }

    cost_pelt_xptr <- fastcpd_test_make_cost_pelt_xptr()
    attr(cost_pelt_xptr, "fastcpd_cost_arity") <- 1L

    result_r <- fastcpd(
      ~ . - 1, data.frame(x = data), family = "custom", cost = cost_pelt_r
    )
    result_xptr <- fastcpd(
      ~ . - 1, data.frame(x = data), family = "custom", cost = cost_pelt_xptr
    )

    testthat::expect_equal(result_xptr@cp_set, result_r@cp_set)
    testthat::expect_equal(
      unname(unlist(result_xptr@cost_values)),
      unname(unlist(result_r@cost_values))
    )
  }
)

testthat::test_that(
  "compiled `cost_gradient`/`cost_hessian` match their R-closure equivalents", {
    # `family = "custom"` drives its sequential gradient descent through
    # `cost_gradient` / `cost_hessian` (an R-closure `cost` alone, with
    # neither, raises `bad_function_call` -- this is a pre-existing
    # restriction of the SeGD path, not specific to compiled costs), so the
    # SeGD-style compiled-cost equivalence is exercised by pairing a
    # (required) R-closure `cost` with compiled `cost_gradient` /
    # `cost_hessian` -- the one compiled/closure mix `check_cost` allows.
    cost_sen_r <- function(data, theta) {
      sum((data[, 1] - as.vector(theta))^2) / 2
    }
    cost_gradient_r <- function(data, theta) {
      -sum(data[, 1] - as.vector(theta))
    }
    cost_hessian_r <- function(data, theta) {
      matrix(nrow(data), 1, 1)
    }

    cost_gradient_xptr <- fastcpd_test_make_cost_gradient_xptr()
    cost_hessian_xptr <- fastcpd_test_make_cost_hessian_xptr()

    result_r <- fastcpd(
      ~ . - 1, data.frame(x = data), family = "custom", cost = cost_sen_r,
      cost_gradient = cost_gradient_r, cost_hessian = cost_hessian_r,
      p = 1
    )
    result_mixed <- fastcpd(
      ~ . - 1, data.frame(x = data), family = "custom", cost = cost_sen_r,
      cost_gradient = cost_gradient_xptr, cost_hessian = cost_hessian_xptr,
      p = 1
    )

    testthat::expect_equal(result_mixed@cp_set, result_r@cp_set)
    testthat::expect_equal(
      unname(unlist(result_mixed@cost_values)),
      unname(unlist(result_r@cost_values)),
      tolerance = 1e-6
    )
  }
)

testthat::test_that(
  "compiled cost without `fastcpd_cost_arity` attribute errors clearly", {
    cost_pelt_xptr <- fastcpd_test_make_cost_pelt_xptr()

    testthat::expect_error(
      fastcpd(
        ~ . - 1, data.frame(x = data), family = "custom", cost = cost_pelt_xptr
      ),
      "fastcpd_cost_arity"
    )
  }
)

testthat::test_that(
  "compiled `cost` can be combined with R-closure `cost_gradient`/`cost_hessian`", {
    # The BFGS warm-start path uses NumericalGradient(cost_function_sen_) so
    # it only needs the cost itself -- no R closure required. The gradient and
    # Hessian closures drive the SeGD parameter update independently and are
    # wrapped by their own std::function adapters. The combination that was
    # previously rejected (XPtr cost + R closure gradient/Hessian) now works.
    cost_sen_xptr <- fastcpd_test_make_cost_sen_xptr()
    attr(cost_sen_xptr, "fastcpd_cost_arity") <- 2L
    cost_gradient_r <- function(data, theta) -sum(data[, 1] - as.vector(theta))
    cost_hessian_r <- function(data, theta) matrix(nrow(data), 1, 1)
    cost_sen_r <- function(data, theta) sum((data[, 1] - as.vector(theta))^2) / 2

    result_xptr <- fastcpd(
      ~ . - 1, data.frame(x = data), family = "custom", cost = cost_sen_xptr,
      cost_gradient = cost_gradient_r, cost_hessian = cost_hessian_r,
      p = 1
    )
    result_r <- fastcpd(
      ~ . - 1, data.frame(x = data), family = "custom", cost = cost_sen_r,
      cost_gradient = cost_gradient_r, cost_hessian = cost_hessian_r,
      p = 1
    )

    testthat::expect_equal(result_xptr@cp_set, result_r@cp_set)
  }
)
