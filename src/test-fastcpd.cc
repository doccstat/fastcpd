#include "fastcpd_test.h"
#include "fastcpd_test_constants.h"
#include "testthat.h"

using ::arma::approx_equal;
using ::arma::colvec;
using ::arma::mat;
using ::arma::ones;

using ::fastcpd::test::FastcpdTest;

context("GetNllPelt Unit Test") {
  test_that("arma(3, 2) is correct for 200 data points") {
    const std::tuple<mat, mat, double> cost_result = FastcpdTest::GetNllPelt(
        colvec(kARMA32.data(), kARMA32.size()), 0, 199, false, R_NilValue);
    const colvec par = std::get<0>(cost_result);
    const colvec residuals = std::get<1>(cost_result);
    const double value = std::get<2>(cost_result);

    // Expected values obtained from the following R code:
    //   arima(x, order = c(3, 0, 2), include.mean = FALSE)
    const colvec expected_par(kARMA32PAR.data(), kARMA32PAR.size());
    const double expected_value(kARMA32VALUE);
    const colvec expected_residuals(kARMA32RESIDUALS.data(),
                                    kARMA32RESIDUALS.size());

    expect_true(norm(par - expected_par, "fro") < 1e-6);
    expect_true(abs(value - expected_value) < 1e-4);
    expect_true(norm(residuals - expected_residuals, "fro") < 1e-6);
  }
}

context("GetNllSen Unit Test") {
  test_that("arma(3, 2) is correct for 200 data points") {
    const double value = FastcpdTest::GetNllSen(
        colvec(kARMA32.data(), kARMA32.size()), 0, 199, 0.1 * ones<colvec>(6));
    const double expected_value = 1363.288;
    expect_true(abs(value - expected_value) < 0.001);
  }
}

context("GetGradient Unit Test") {
  test_that("arma(3, 2) is correct for 200 data points") {
    const colvec gradient = FastcpdTest::GetGradientArma(
        colvec(kARMA32.data(), kARMA32.size()), 0, 199, 0.1 * ones<colvec>(6));
    const colvec expected_gradient = {4.401258, 6.600128, -7.591818,
                                      4.151778, 7.503752, -2.806806};
    expect_true(norm(gradient - expected_gradient, "fro") < 1e-6);
  }
}

context("GetHessian Unit Test") {
  test_that("binomial is correct for a two-dimensional data") {
    const mat hessian =
        FastcpdTest::GetHessianBinomial({{1, 1, 0.2}}, 0, 0, {-0.5, 0.3});
    const mat expected_hessian = {{0.23828, 0.047656}, {0.047656, 0.0095312}};
    expect_true(norm(hessian - expected_hessian, "fro") < 1e-6);
  }

  test_that("poisson is correct for a two-dimensional data") {
    const mat hessian =
        FastcpdTest::GetHessianPoisson({{4, 1, 0.2}}, 0, 0, {-0.5, 0.3});
    const mat expected_hessian = {{0.6440364, 0.128807}, {0.128807, 0.0257616}};
    expect_true(norm(hessian - expected_hessian, "fro") < 1e-6);
  }

  test_that("arma(3, 2) is correct for 200 data points") {
    const mat hessian = FastcpdTest::GetHessianArma(
        colvec(kARMA32.data(), kARMA32.size()), 0, 199, 0.1 * ones<colvec>(6));
    const mat expected_hessian = {
        {12.406525, 18.60483, -21.40027, 4.743794, 28.98263, -44.01258},
        {18.604831, 27.89981, -32.09185, 25.380851, 27.48253, -66.00128},
        {-21.400265, -32.09185, 36.91375, -24.424268, -34.63643, 75.91818},
        {4.743794, 25.38085, -24.42427, -4.424596, 35.77934, -41.51778},
        {28.982631, 27.48253, -34.63643, 35.779338, 24.80587, -75.03752},
        {-44.012575, -66.00128, 75.91818, -41.51778, -75.03752, 106.13612}};
    expect_true(norm(hessian - expected_hessian, "fro") < 2e-5);
  }
}
