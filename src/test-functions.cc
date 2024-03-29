#include "fastcpd_classes.h"
#include "testthat.h"

// ``` r
// set.seed(1)
// x <- arima.sim(list(ar = c(0.1, -0.2, 0.6), ma = c(0.1, 0.5)), n = 200)
// ```
constexpr std::array<double, 200> kARMA32 = {
   0.16086394,  1.36633450,  0.71425776, -0.23698470,  1.13222067,
  -0.81209066,  1.04080322,  2.50485157,  0.10240048,  0.04309484,
   1.76843576, -0.37048506,  2.30809336,  1.49936227,  1.35265076,
   1.29760473,  0.36324437,  0.71686001, -1.38011604,  1.31601817,
   0.23507054,  1.85294823,  1.69728898,  0.36409193,  1.58620681,
  -0.12382253, -1.15284844,  0.56023783, -0.82866400, -0.78412375,
   0.27626904, -0.89428088, -1.20560631, -0.26274899,  0.55851260,
  -2.08831014,  0.55245008,  0.43857619,  0.07374949,  0.21973100,
   1.14151925,  0.26646413, -0.40061958,  1.87872196,  1.43780397,
   0.94785207,  3.17151372,  2.05753994, -0.28716435,  1.04100239,
  -0.52417819, -1.31541151, -0.68211713, -0.37625744, -1.90734262,
  -0.43675851, -0.98220438,  0.62556888,  0.56308654,  1.20736086,
   1.21701663,  2.39376315,  0.44496334,  0.61182334,  2.47669709,
  -0.34600745,  0.28084458,  0.84441002, -0.64229642, -0.57212145,
   0.88417372, -0.45000449, -0.84170441,  1.74011253, -0.26090482,
  -0.40863069,  0.82994931,  0.62104715, -0.40145671,  0.64495703,
  -0.20479077, -0.80159791,  0.03467575, -0.70410088, -0.05561702,
  -1.60362778, -0.15123335, -1.99275907, -1.43254808, -1.16204543,
  -1.88076120, -1.20144209, -2.68731572, -0.20020090, -2.70791602,
  -1.88487685, -1.76178438, -2.51812889,  0.42374247, -0.66036783,
  -1.90264031, -1.56448276, -0.52916463, -1.67544799, -1.09523234,
  -1.06238132, -2.63189117, -2.39647166, -0.20191830, -2.17889216,
  -2.56192766,  1.47864928, -0.72726882, -1.16710509,  2.16310135,
   0.88631727, -1.04561085,  3.60342576,  0.75721680, -1.61927610,
   1.43432190,  0.40246854,  1.03834924,  1.32430373,  1.78780989,
   0.55837792,  0.37964132,  0.89227190,  0.96549226,  2.28252559,
   2.19925399,  2.69098578,  0.60192677,  2.30626534,  1.42748530,
  -0.91065824,  1.49145646,  0.34749113,  0.89928610,  0.21630369,
   0.27574153, -0.82828660, -0.49301554,  0.20304732, -1.15816310,
   0.50596914, -1.08686163, -1.65838396,  1.08743329, -1.60742076,
  -0.34229029,  0.09191278, -0.30952153,  1.28456656,  2.20568531,
   0.45012999, -1.15041623,  3.22202770,  0.59657857,  0.58186733,
   2.24635394,  1.24564121, -0.09565102,  1.74843029,  0.50099276,
  -1.55686262,  1.44386747,  1.68977984, -0.71676002, -0.06760279,
   1.51300168, -0.87236517, -1.84696719,  0.70403277, -1.58031874,
  -1.80995143, -0.94833112,  0.83596631, -1.54181203, -1.62742880,
   0.51827539, -1.06763563, -2.04778834, -2.53620199, -1.90693103,
  -1.85012658, -1.64826101, -0.75785666, -0.33506819, -2.05035803,
   0.66580539, -0.21328442, -0.12953955,  1.53135295,  0.73161908,
  -0.65013812,  1.89430814, -1.56479241, -1.08874870,  0.03043624
};

// ``` r
// arima(x, order = c(3, 0, 2), include.mean = FALSE)
// ```
constexpr std::array<double, 6> kARMA32PAR =
  {0.09403583, -0.21094577, 0.64057694, 0.08878535, 0.50335714, 0.976383};
constexpr double kARMA32VALUE = 282.2705;
constexpr std::array<double, 200> kARMA32RESIDUALS = {
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

context("negative_log_likelihood_wo_theta Unit Test") {
  test_that("arma(3, 2) is correct for 200 data points") {
    fastcpd::classes::Fastcpd fastcpd_class(
      /* beta */ 0,
      /* convexity_coef */ 0,
      /* cost */ R_NilValue,
      /* cost_adjustment */ "MBIC",
      /* cost_gradient */ R_NilValue,
      /* cost_hessian */ R_NilValue,
      /* cp_only */ true,
      /* data */ mat(),
      /* epsilon */ 0.0,
      /* family */ "arma",
      /* k */ R_NilValue,
      /* line_search */ colvec(),
      /* lower */ colvec(),
      /* momentum_coef */ 0.0,
      /* order */ colvec({3, 2}),
      /* p */ 0,
      /* p_response */ 0,
      /* r_progress */ false,
      /* segment_count */ 0,
      /* trim */ 0,
      /* upper */ colvec(),
      /* vanilla_percentage */ 0.0,
      /* variance_estimate */ mat(),
      /* warm_start */ false
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
      /* beta */ 0,
      /* convexity_coef */ 0,
      /* cost */ R_NilValue,
      /* cost_adjustment */ "MBIC",
      /* cost_gradient */ R_NilValue,
      /* cost_hessian */ R_NilValue,
      /* cp_only */ true,
      /* data */ mat(),
      /* epsilon */ 0.0,
      /* family */ "arma",
      /* k */ R_NilValue,
      /* line_search */ colvec(),
      /* lower */ colvec(),
      /* momentum_coef */ 0.0,
      /* order */ colvec({3, 2}),
      /* p */ 0,
      /* p_response */ 0,
      /* r_progress */ false,
      /* segment_count */ 0,
      /* trim */ 0,
      /* upper */ colvec(),
      /* vanilla_percentage */ 0.0,
      /* variance_estimate */ mat(),
      /* warm_start */ false
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
      /* beta */ 0,
      /* convexity_coef */ 0,
      /* cost */ R_NilValue,
      /* cost_adjustment */ "MBIC",
      /* cost_gradient */ R_NilValue,
      /* cost_hessian */ R_NilValue,
      /* cp_only */ true,
      /* data */ mat(),
      /* epsilon */ 0.0,
      /* family */ "arma",
      /* k */ R_NilValue,
      /* line_search */ colvec(),
      /* lower */ colvec(),
      /* momentum_coef */ 0.0,
      /* order */ colvec({3, 2}),
      /* p */ 0,
      /* p_response */ 0,
      /* r_progress */ false,
      /* segment_count */ 0,
      /* trim */ 0,
      /* upper */ colvec(),
      /* vanilla_percentage */ 0.0,
      /* variance_estimate */ mat(),
      /* warm_start */ false
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
      /* beta */ 0,
      /* convexity_coef */ 0,
      /* cost */ R_NilValue,
      /* cost_adjustment */ "MBIC",
      /* cost_gradient */ R_NilValue,
      /* cost_hessian */ R_NilValue,
      /* cp_only */ true,
      /* data */ mat(),
      /* epsilon */ 0.0,
      /* family */ "binomial",
      /* k */ R_NilValue,
      /* line_search */ colvec(),
      /* lower */ colvec(),
      /* momentum_coef */ 0.0,
      /* order */ colvec(),
      /* p */ 0,
      /* p_response */ 0,
      /* r_progress */ false,
      /* segment_count */ 0,
      /* trim */ 0,
      /* upper */ colvec(),
      /* vanilla_percentage */ 0.0,
      /* variance_estimate */ mat(),
      /* warm_start */ false
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
      /* beta */ 0,
      /* convexity_coef */ 0,
      /* cost */ R_NilValue,
      /* cost_adjustment */ "MBIC",
      /* cost_gradient */ R_NilValue,
      /* cost_hessian */ R_NilValue,
      /* cp_only */ true,
      /* data */ mat(),
      /* epsilon */ 0.0,
      /* family */ "poisson",
      /* k */ R_NilValue,
      /* line_search */ colvec(),
      /* lower */ colvec(),
      /* momentum_coef */ 0.0,
      /* order */ colvec(),
      /* p */ 0,
      /* p_response */ 0,
      /* r_progress */ false,
      /* segment_count */ 0,
      /* trim */ 0,
      /* upper */ colvec(),
      /* vanilla_percentage */ 0.0,
      /* variance_estimate */ mat(),
      /* warm_start */ false
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
      /* beta */ 0,
      /* convexity_coef */ 0,
      /* cost */ R_NilValue,
      /* cost_adjustment */ "MBIC",
      /* cost_gradient */ R_NilValue,
      /* cost_hessian */ R_NilValue,
      /* cp_only */ true,
      /* data */ mat(),
      /* epsilon */ 0.0,
      /* family */ "arma",
      /* k */ R_NilValue,
      /* line_search */ colvec(),
      /* lower */ colvec(),
      /* momentum_coef */ 0.0,
      /* order */ colvec({3, 2}),
      /* p */ 0,
      /* p_response */ 0,
      /* r_progress */ false,
      /* segment_count */ 0,
      /* trim */ 0,
      /* upper */ colvec(),
      /* vanilla_percentage */ 0.0,
      /* variance_estimate */ mat(),
      /* warm_start */ false
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

context("update_theta_sum Unit Test") {
  test_that("update performs normally") {
    fastcpd::classes::Fastcpd fastcpd_class(
      /* beta */ 0,
      /* convexity_coef */ 0,
      /* cost */ R_NilValue,
      /* cost_adjustment */ "MBIC",
      /* cost_gradient */ R_NilValue,
      /* cost_hessian */ R_NilValue,
      /* cp_only */ true,
      /* data */ mat(),
      /* epsilon */ 0.0,
      /* family */ "gaussian",
      /* k */ R_NilValue,
      /* line_search */ colvec(),
      /* lower */ colvec(),
      /* momentum_coef */ 0.0,
      /* order */ colvec(),
      /* p */ 3,
      /* p_response */ 0,
      /* r_progress */ false,
      /* segment_count */ 0,
      /* trim */ 0,
      /* upper */ colvec(),
      /* vanilla_percentage */ 0.0,
      /* variance_estimate */ mat(),
      /* warm_start */ false
    );
    fastcpd_class.create_theta_sum(0, colvec({1, 2, 3}));
    fastcpd_class.update_theta_sum(0, colvec({4, 5, 6}));
    expect_true(fastcpd_class.get_theta_sum().n_rows == 3);
    expect_true(fastcpd_class.get_theta_sum().n_cols == 1);
    colvec theta_sum = fastcpd_class.get_theta_sum().col(0);
    colvec expected_theta_sum = {5, 7, 9};
    expect_true(approx_equal(theta_sum, expected_theta_sum, "absdiff", 1e-6));
  }
}
