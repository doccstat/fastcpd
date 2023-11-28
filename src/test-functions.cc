#include "fastcpd_classes.h"
#include "fastcpd_constants.h"
#include "testthat.h"

context("negative_log_likelihood_wo_theta Unit Test") {
  test_that("arma(3, 2) is correct for 200 data points") {
    fastcpd::classes::Fastcpd fastcpd_class(
      /* data */ mat(),
      /* beta */ 0,
      /* p */ 0,
      /* order */ colvec({3, 2}),
      /* family */ "arma",
      /* vanilla_percentage */ 0.0,
      /* segment_count */ 0,
      /* winsorise_minval */ 0.0,
      /* winsorise_maxval */ 0.0,
      /* epsilon */ 0.0,
      /* min_prob */ 0.0,
      /* momentum_coef */ 0.0,
      /* lower */ colvec(),
      /* upper */ colvec(),
      /* mean_data_cov */ mat()
    );

    const colvec data(kARMA32.data(), kARMA32.size());
    const List out = fastcpd_class.negative_log_likelihood_wo_theta(
      data, 0, false, R_NilValue
    );
    const colvec par = out["par"];
    const double value = out["value"];
    const colvec residuals = out["residuals"];

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

context("negative_log_likelihood_wo_cv Unit Test") {
  test_that("arma(3, 2) is correct for 200 data points") {
    fastcpd::classes::Fastcpd fastcpd_class(
      /* data */ mat(),
      /* beta */ 0,
      /* p */ 0,
      /* order */ colvec({3, 2}),
      /* family */ "arma",
      /* vanilla_percentage */ 0.0,
      /* segment_count */ 0,
      /* winsorise_minval */ 0.0,
      /* winsorise_maxval */ 0.0,
      /* epsilon */ 0.0,
      /* min_prob */ 0.0,
      /* momentum_coef */ 0.0,
      /* lower */ colvec(),
      /* upper */ colvec(),
      /* mean_data_cov */ mat()
    );

    const colvec data(kARMA32.data(), kARMA32.size());
    const colvec theta = 0.1 * ones<colvec>(6);
    const double value =
      fastcpd_class.negative_log_likelihood_wo_cv(data, theta, 0.0);
    const double expected_value = 1363.288;
    expect_true(abs(value -  expected_value) < 0.001);
  }
}

context("cost_update_gradient Unit Test") {
  test_that("arma(3, 2) is correct for 200 data points") {
    fastcpd::classes::Fastcpd fastcpd_class(
      /* data */ mat(),
      /* beta */ 0,
      /* p */ 0,
      /* order */ colvec({3, 2}),
      /* family */ "arma",
      /* vanilla_percentage */ 0.0,
      /* segment_count */ 0,
      /* winsorise_minval */ 0.0,
      /* winsorise_maxval */ 0.0,
      /* epsilon */ 0.0,
      /* min_prob */ 0.0,
      /* momentum_coef */ 0.0,
      /* lower */ colvec(),
      /* upper */ colvec(),
      /* mean_data_cov */ mat()
    );

    const colvec data(kARMA32.data(), kARMA32.size());
    const colvec theta = 0.1 * ones<colvec>(6);
    const colvec gradient = fastcpd_class.cost_update_gradient(data, theta);
    const colvec expected_gradient =
      {4.401258, 6.600128, -7.591818, 4.151778, 7.503752, -2.806806};
    expect_true(norm(gradient - expected_gradient, "fro") < 1e-6);
  }
}

context("cost_update_hessian Unit Test") {
  test_that("binomal is correct for a two dimensional data") {
    fastcpd::classes::Fastcpd fastcpd_class(
      /* data */ mat(),
      /* beta */ 0,
      /* p */ 0,
      /* order */ colvec(),
      /* family */ "binomial",
      /* vanilla_percentage */ 0.0,
      /* segment_count */ 0,
      /* winsorise_minval */ 0.0,
      /* winsorise_maxval */ 0.0,
      /* epsilon */ 0.0,
      /* min_prob */ 0.0,
      /* momentum_coef */ 0.0,
      /* lower */ colvec(),
      /* upper */ colvec(),
      /* mean_data_cov */ mat()
    );

    const mat data = {{1, 1, 0.2}};
    const colvec theta = {-0.5, 0.3};
    const mat hessian = fastcpd_class.cost_update_hessian(data, theta);
    const mat expected_hessian =
        {{0.238'28, 0.047'656}, {0.047'656, 0.009'531'2}};
    expect_true(norm(hessian - expected_hessian, "fro") < 0.000'001);
  }

  test_that("poisson is correct for a two dimensional data") {
    fastcpd::classes::Fastcpd fastcpd_class(
      /* data */ mat(),
      /* beta */ 0,
      /* p */ 0,
      /* order */ colvec(),
      /* family */ "poisson",
      /* vanilla_percentage */ 0.0,
      /* segment_count */ 0,
      /* winsorise_minval */ 0.0,
      /* winsorise_maxval */ 0.0,
      /* epsilon */ 0.0,
      /* min_prob */ 1e10,
      /* momentum_coef */ 0.0,
      /* lower */ colvec(),
      /* upper */ colvec(),
      /* mean_data_cov */ mat()
    );

    const mat data = {{4, 1, 0.2}};
    const colvec theta = {-0.5, 0.3};
    const mat hessian = fastcpd_class.cost_update_hessian(data, theta);
    const mat expected_hessian =
        {{0.644'036'4, 0.128'807}, {0.128'807, 0.025'761'6}};
    expect_true(norm(hessian - expected_hessian, "fro") < 0.000'001);
  }

  test_that("arma(3, 2) is correct for 200 data points") {
    fastcpd::classes::Fastcpd fastcpd_class(
      /* data */ mat(),
      /* beta */ 0,
      /* p */ 0,
      /* order */ colvec({3, 2}),
      /* family */ "arma",
      /* vanilla_percentage */ 0.0,
      /* segment_count */ 0,
      /* winsorise_minval */ 0.0,
      /* winsorise_maxval */ 0.0,
      /* epsilon */ 0.0,
      /* min_prob */ 0.0,
      /* momentum_coef */ 0.0,
      /* lower */ colvec(),
      /* upper */ colvec(),
      /* mean_data_cov */ mat()
    );

    const colvec data(kARMA32.data(), kARMA32.size());
    const colvec theta = 0.1 * ones<colvec>(6);
    const mat hessian = fastcpd_class.cost_update_hessian(data, theta);
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
