#include "fastcpd.h"

//' Function to calculate the Hessian matrix at the current data.
//'
//' @param data A data frame containing the data to be segmented.
//' @param theta Estimated theta from the previous iteration.
//' @param family Family of the model.
//' @param min_prob Minimum probability to avoid numerical issues.
//'
//' @return Hessian at the current data.
// [[Rcpp::export]]
arma::mat cost_update_hessian(arma::mat data,
                              arma::colvec theta,
                              std::string family,
                              double min_prob) {
    arma::rowvec new_data = data.row(data.n_rows - 1);
    arma::rowvec x = new_data.tail(new_data.n_elem - 1);
    arma::mat hessian;
    if (family.compare("binomial") == 0) {
        double prob = 1 / (1 + exp(-arma::as_scalar(x * theta)));
        hessian = (x.t() * x) * arma::as_scalar((1 - prob) * prob);
    } else if (family.compare("poisson") == 0) {
        double prob = exp(arma::as_scalar(x * theta));
        hessian = (x.t() * x) * std::min(arma::as_scalar(prob), min_prob);
    } else {
        // `family` is either "lasso" or "gaussian".
        hessian = x.t() * x;
    }
    return hessian;
}

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
}
