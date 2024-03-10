#include "fastcpd_classes.h"

namespace fastcpd::classes {

CostFunction::CostFunction(Function cost) : cost(cost) {}

List CostFunction::operator()(  // # nocov
    mat data,
    Nullable<colvec> theta,
    double lambda,  // UNUSED
    bool cv,  // UNUSED
    Nullable<colvec> start  // UNUSED
) {
  DEBUG_RCOUT(data.n_rows);
  SEXP value =
    theta.isNull() ? cost(data) : cost(data, as<colvec>(theta));  // # nocov
  return List::create(
      Named("value") = as<double>(value)  // # nocov
  );
}

CostGradient::CostGradient(Function cost_gradient) :
    cost_gradient(cost_gradient) {}

colvec CostGradient::operator()(mat data, colvec theta) {
  return as<colvec>(cost_gradient(data, theta));
}

CostHessian::CostHessian(Function cost_hessian) :
    cost_hessian(cost_hessian) {}

mat CostHessian::operator()(mat data, colvec theta) {
  return as<mat>(cost_hessian(data, theta));
}

}  // namespace fastcpd::classes
