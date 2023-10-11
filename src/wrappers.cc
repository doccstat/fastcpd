#include "wrappers.h"

using ::Rcpp::as;

namespace fastcpd::wrappers {

CostFunction::CostFunction(Rcpp::Function cost) : cost(cost) {}

List CostFunction::operator()(
    arma::mat data,
    Nullable<arma::colvec> theta,
    std::string family,  // UNUSED
    double lambda,  // UNUSED
    bool cv,  // UNUSED
    Nullable<arma::colvec> start  // UNUSED
) {
  return theta.isNull()? cost(data) : cost(data, theta);
}

CostGradient::CostGradient(Rcpp::Function cost_gradient) :
    cost_gradient(cost_gradient) {}

arma::colvec CostGradient::operator()(
    arma::mat data,
    arma::colvec theta,
    std::string family  // UNUSED
) {
  return as<arma::colvec>(cost_gradient(data, theta));
}

CostHessian::CostHessian(Rcpp::Function cost_hessian) :
    cost_hessian(cost_hessian) {}

arma::mat CostHessian::operator()(
    arma::mat data,
    arma::colvec theta,
    std::string family,  // UNUSED
    double min_prob  // UNUSED
) {
  return as<arma::mat>(cost_hessian(data, theta));
}

}  // namespace fastcpd::wrappers
