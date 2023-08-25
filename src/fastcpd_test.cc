#include "fastcpd.h"

context("cost_update_hessian Unit Test") {
  test_that("binomal is correct for a two dimensional data") {
    arma::colvec theta = {-0.5, 0.3};
    arma::mat data = {{1, 1, 0.2}};
    std::string family = "binomial";
    double min_prob = 0.0;
    arma::mat hessian = cost_update_hessian(data, theta, family, min_prob);
    arma::mat hessian_expected = {{0.238280, 0.0476560}, {0.047656, 0.0095312}};
    expect_true(arma::norm(hessian - hessian_expected, "fro") < 0.000001);
  }

  test_that("poisson is correct for a two dimensional data") {
    arma::colvec theta = {-0.5, 0.3};
    arma::mat data = {{4, 1, 0.2}};
    std::string family = "poisson";
    double min_prob = 1e10;
    arma::mat hessian = cost_update_hessian(data, theta, family, min_prob);
    arma::mat hessian_expected = {{0.6440364, 0.128807}, {0.128807, 0.0257616}};
    expect_true(arma::norm(hessian - hessian_expected, "fro") < 0.000001);
  }
}
