// Commented out due to the LTO on CRAN.
// #include "fastcpd.h"
// #include "testthat.h"

// context("cost_update_hessian Unit Test") {
//   test_that("binomal is correct for a two dimensional data") {
//     arma::colvec theta = {-0.5, 0.3};
//     arma::mat data = {{1, 1, 0.2}};
//     std::string family = "binomial";
//     double min_prob = 0.0;
//     arma::mat hessian = cost_update_hessian(data, theta, family, min_prob);
//     arma::mat hessian_expected =
//         {{0.238'28, 0.047'656}, {0.047'656, 0.009'531'2}};
//     expect_true(arma::norm(hessian - hessian_expected, "fro") < 0.000'001);
//   }

//   test_that("poisson is correct for a two dimensional data") {
//     arma::colvec theta = {-0.5, 0.3};
//     arma::mat data = {{4, 1, 0.2}};
//     std::string family = "poisson";
//     double min_prob = 1e10;
//     arma::mat hessian = cost_update_hessian(data, theta, family, min_prob);
//     arma::mat hessian_expected =
//         {{0.644'036'4, 0.128'807}, {0.128'807, 0.025'761'6}};
//     expect_true(arma::norm(hessian - hessian_expected, "fro") < 0.000'001);
//   }
// }
