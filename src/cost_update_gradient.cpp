#include "fastcpd.h"

//' Function to calculate the gradient at the current data.
//'
//' @param data A data frame containing the data to be segmented.
//' @param theta Estimated theta from the previous iteration.
//' @param family Family of the model.
//'
//' @return Gradient at the current data.
// [[Rcpp::export]]
arma::colvec cost_update_gradient(arma::mat data,
                                  arma::colvec theta,
                                  std::string family) {
    arma::rowvec new_data = data.row(data.n_rows - 1);
    arma::rowvec x = new_data.tail(new_data.n_elem - 1);
    double y = new_data(0);
    arma::colvec gradient;
    if (family.compare("binomial") == 0) {
        gradient = - (y - 1 / (1 + exp(-arma::as_scalar(x * theta)))) * x.t();
    } else if (family.compare("poisson") == 0) {
        gradient = - (y - exp(arma::as_scalar(x * theta))) * x.t();
    } else {
        // `family` is either "lasso" or "gaussian".
        gradient = - (y - arma::as_scalar(x * theta)) * x.t();
    }
    return gradient;
}
