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
