#include "fastcpd_test.h"
#include "fastcpd_test_constants.h"
#include "testthat.h"

using ::arma::approx_equal;
using ::arma::ones;

using ::fastcpd::classes::CostResult;
using ::fastcpd::test::FastcpdTest;

context("get_nll_pelt Unit Test") {
  test_that("arma(3, 2) is correct for 200 data points") {
    const CostResult cost_result = FastcpdTest::get_nll_pelt(
      colvec(kARMA32.data(), kARMA32.size()), 0, 199, 0, false, R_NilValue
    );
    const colvec par = cost_result.par;
    const double value = cost_result.value;
    const colvec residuals = cost_result.residuals;

    // Expected values obtained from the following R code
    // arima(x, order = c(3, 0, 2), include.mean = FALSE)
    const colvec expected_par(kARMA32PAR.data(), kARMA32PAR.size());
    const double expected_value(kARMA32VALUE);
    const colvec expected_residuals(
      kARMA32RESIDUALS.data(), kARMA32RESIDUALS.size()
    );

    expect_true(norm(par - expected_par, "fro") < 0.000'001);
    expect_true(abs(value - expected_value) < 1e-4);
    expect_true(norm(residuals - expected_residuals, "fro") < 0.000'001);
  }
}

context("get_nll_sen Unit Test") {
  test_that("arma(3, 2) is correct for 200 data points") {
    const double value = FastcpdTest::get_nll_sen(
      colvec(kARMA32.data(), kARMA32.size()), 0, 199, 0.1 * ones<colvec>(6), 0.0
    );
    const double expected_value = 1363.288;
    expect_true(abs(value -  expected_value) < 0.001);
  }
}

context("get_gradient Unit Test") {
  test_that("arma(3, 2) is correct for 200 data points") {
    const colvec gradient = FastcpdTest::get_gradient_arma(
      colvec(kARMA32.data(), kARMA32.size()), 0, 199, 0.1 * ones<colvec>(6)
    );
    const colvec expected_gradient =
      {4.401258, 6.600128, -7.591818, 4.151778, 7.503752, -2.806806};
    expect_true(norm(gradient - expected_gradient, "fro") < 1e-6);
  }
}

context("get_hessian Unit Test") {
  test_that("binomal is correct for a two dimensional data") {
    const mat hessian = FastcpdTest::get_hessian_binomial(
      {{1, 1, 0.2}}, 0, 0, {-0.5, 0.3}
    );
    const mat expected_hessian =
        {{0.238'28, 0.047'656}, {0.047'656, 0.009'531'2}};
    expect_true(norm(hessian - expected_hessian, "fro") < 0.000'001);
  }

  test_that("poisson is correct for a two dimensional data") {
    const mat hessian = FastcpdTest::get_hessian_poisson(
      {{4, 1, 0.2}}, 0, 0, {-0.5, 0.3}
    );
    const mat expected_hessian =
        {{0.644'036'4, 0.128'807}, {0.128'807, 0.025'761'6}};
    expect_true(norm(hessian - expected_hessian, "fro") < 0.000'001);
  }

  test_that("arma(3, 2) is correct for 200 data points") {
    const mat hessian = FastcpdTest::get_hessian_arma(
      colvec(kARMA32.data(), kARMA32.size()), 0, 199, 0.1 * ones<colvec>(6)
    );
    const mat expected_hessian = {
      { 12.406525,  18.60483, -21.40027,   4.743794,  28.98263, -44.01258},
      { 18.604831,  27.89981, -32.09185,  25.380851,  27.48253, -66.00128},
      {-21.400265, -32.09185,  36.91375, -24.424268, -34.63643,  75.91818},
      {  4.743794,  25.38085, -24.42427,  -4.424596,  35.77934, -41.51778},
      { 28.982631,  27.48253, -34.63643,  35.779338,  24.80587, -75.03752},
      {-44.012575, -66.00128,  75.91818, -41.517780, -75.03752, 106.13612}
    };
    expect_true(norm(hessian - expected_hessian, "fro") < 2e-5);
  }
}

context("update_theta_sum Unit Test") {
  test_that("update performs normally") {
    mat theta_sum = FastcpdTest::update_theta_sum(0, {1, 2, 3}, {4, 5, 6});
    expect_true(theta_sum.n_rows == 3);
    expect_true(theta_sum.n_cols == 1);
    colvec theta_sum_col = theta_sum.col(0);
    colvec expected_theta_sum = {5, 7, 9};
    expect_true(
      approx_equal(theta_sum_col, expected_theta_sum, "absdiff", 1e-6)
    );
  }
}
