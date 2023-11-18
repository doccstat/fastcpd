#include "constants.h"
#include "functions.h"
#include "testthat.h"

using ::fastcpd::functions::negative_log_likelihood_wo_theta;
using ::fastcpd::functions::negative_log_likelihood_wo_cv;
using ::fastcpd::functions::cost_update_gradient;
using ::fastcpd::functions::cost_update_hessian;

context("negative_log_likelihood_wo_theta Unit Test") {
  test_that("arma(3, 2) is correct for 200 data points") {
    colvec data(time_series_arma_3_2.data(), time_series_arma_3_2.size());
    const colvec order = {3, 2};
    List out = negative_log_likelihood_wo_theta(
      data, "arma", 0, false, R_NilValue, order
    );
    colvec par = out["par"];
    double value = out["value"];
    colvec residuals = out["residuals"];

    // Expected values obtained from the following R code
    // arima(x, order = c(3, 0, 2), include.mean = FALSE)
    colvec expected_par =
      {0.09403583, -0.21094577, 0.64057694, 0.08878535, 0.50335714, 0.976383};
    double expected_value = 282.2705;
    colvec expected_residuals = {
       0.119657519,  1.032474096,  0.272378503, -0.513466262,  0.360896651,
      -1.192879019,  1.433102017,  1.973919976, -0.286092253, -1.068902404,
       0.420311160, -0.092947649,  2.484881682, -0.102348441,  0.694299134,
      -0.001940191, -0.783191507,  0.160458807, -1.822138707,  1.445339739,
       0.149843828,  2.251696352,  0.454279024, -0.728963510,  0.559105842,
      -0.966134762, -1.235481261,  0.222442991, -0.443078439,  0.077838894,
       0.032441912, -0.596905740, -0.524275873, -0.167990371,  1.180571667,
      -1.444230828,  0.568930225,  0.264786828,  1.176883663, -0.286347443,
       0.288502676,  0.276749743, -0.495423915,  1.146055554,  1.153559236,
       0.786288702,  1.531751169,  0.506444241, -1.234784305, -0.674854658,
      -1.319203758, -0.405755195, -0.635780870,  0.006869776, -0.853813481,
       0.172526800, -0.688002195,  1.821840545,  0.761404294,  0.930911151,
       0.355627614,  1.673152098, -0.624378265, -0.491415139,  1.337553158,
      -0.606275222, -0.175528042, -0.520746009, -0.306225934, -0.224190818,
       0.435620199, -0.168223465, -0.450724104,  1.282649234, -0.200834487,
      -0.105649963, -0.190865851,  0.694058623,  0.011427761, -0.068304761,
      -0.747641695, -0.288364059,  0.055643618, -0.605061321,  0.557105497,
      -1.514038247,  0.292865107, -1.545089404, -0.260046835, -0.550001959,
      -0.617433203, -0.020388020, -1.914096665,  1.184042387, -1.628003140,
      -0.402488130, -1.172315731, -0.708756546,  2.149325178,  0.063084092,
      -1.225579201, -1.719248215,  0.409163368,  0.092147088, -0.261270269,
      -0.997034168, -1.469736265, -1.039148546,  0.980849961, -0.543523639,
      -1.309963958,  1.779170492,  0.490420881, -0.084789517,  0.932919342,
       0.862433833, -0.471201838,  2.110806263, -0.320181777, -1.294621935,
      -0.285838658,  0.117991304,  2.473740888,  0.113743496,  0.369227313,
      -0.085562152, -0.322310247, -0.099185568,  0.775029832,  2.117880710,
       1.038558938,  1.188934964, -1.277658296,  0.923501768,  0.174929192,
      -1.424358847,  0.439282889, -0.221313577,  1.563105813, -0.777733380,
      -0.495242117, -0.949201880, -0.161960697,  0.390218501, -0.703796677,
       0.839591893, -1.229100273, -1.021043877,  1.399329898, -0.973579298,
       0.482653713, -0.464357456,  0.437588307,  1.747210767,  1.585330654,
      -0.308262364, -2.320946623,  2.373482023,  0.720108822,  0.683726565,
      -0.169650185,  0.445895234, -0.065852182,  0.362627283, -0.480577366,
      -1.313740882,  0.934489498,  1.482978448, -0.175840543, -1.299513687,
       0.489615625,  0.040886682, -1.652549176, -0.149360323, -0.632232516,
      -0.198392944, -1.226626417,  1.764426514, -0.200280156, -1.568975360,
       0.050687002,  0.313232606, -0.848895099, -2.983143430, -0.724349868,
       0.671856365,  0.053046331, -0.114494834,  0.557116433, -1.114710260,
       1.091938186, -0.029622117,  0.797372863,  1.016158760,  0.205332961,
      -0.842645368,  1.100285732, -2.022265133, -0.499831533, -0.348419523
    };

    expect_true(norm(par - expected_par, "fro") < 0.000'001);
    expect_true(abs(value - expected_value) < 1e-4);
    expect_true(norm(residuals - expected_residuals, "fro") < 0.000'001);
  }
}

context("negative_log_likelihood_wo_cv Unit Test") {
  test_that("arma(3, 2) is correct for 200 data points") {
    colvec data(time_series_arma_3_2.data(), time_series_arma_3_2.size()),
           theta = 0.1 * ones<colvec>(6);
    const colvec order = {3, 2};
    double value = negative_log_likelihood_wo_cv(
      data, theta, "arma", 0.0, R_NilValue, order
    );
    expect_true(abs(value -  1363.288) < 0.001);
  }
}

context("cost_update_gradient Unit Test") {
  test_that("arma(3, 2) is correct for 200 data points") {
    colvec theta = 0.1 * ones<colvec>(6),
           data(time_series_arma_3_2.data(), time_series_arma_3_2.size());
    const colvec order = {3, 2};
    colvec gradient = cost_update_gradient(data, theta, "arma", order),
           expected_gradient =
      {4.401258, 6.600128, -7.591818, 4.151778, 7.503752, -2.806806};
    expect_true(norm(gradient - expected_gradient, "fro") < 1e-6);
  }
}

context("cost_update_hessian Unit Test") {
  test_that("binomal is correct for a two dimensional data") {
    colvec theta = {-0.5, 0.3};
    mat data = {{1, 1, 0.2}};
    std::string family = "binomial";
    double min_prob = 1e10;
    const colvec order = {3, 2};
    mat hessian = cost_update_hessian(data, theta, family, min_prob, order);
    mat expected_hessian =
        {{0.238'28, 0.047'656}, {0.047'656, 0.009'531'2}};
    expect_true(norm(hessian - expected_hessian, "fro") < 0.000'001);
  }

  test_that("poisson is correct for a two dimensional data") {
    colvec theta = {-0.5, 0.3};
    mat data = {{4, 1, 0.2}};
    std::string family = "poisson";
    double min_prob = 1e10;
    const colvec order = {3, 2};
    mat hessian = cost_update_hessian(data, theta, family, min_prob, order);
    mat expected_hessian =
        {{0.644'036'4, 0.128'807}, {0.128'807, 0.025'761'6}};
    expect_true(norm(hessian - expected_hessian, "fro") < 0.000'001);
  }

  test_that("arma(3, 2) is correct for 200 data points") {
    colvec theta = 0.1 * ones<colvec>(6);
    colvec data(time_series_arma_3_2.data(), time_series_arma_3_2.size());
    const double min_prob = 1e10;
    const colvec order = {3, 2};
    mat hessian = cost_update_hessian(data, theta, "arma", min_prob, order);
    mat expected_hessian = {
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
