#include "wrappers.h"

namespace fastcpd::wrappers {

CostFunction::CostFunction(Function cost) : cost(cost) {}

List CostFunction::operator()(  // # nocov
    mat data,
    Nullable<colvec> theta,
    string family,  // UNUSED
    double lambda,  // UNUSED
    bool cv,  // UNUSED
    Nullable<colvec> start,  // UNUSED
    const colvec order  // UNUSED
) {
  return theta.isNull()? cost(data) : cost(data, theta);  // # nocov
}

CostGradient::CostGradient(Function cost_gradient) :
    cost_gradient(cost_gradient) {}

colvec CostGradient::operator()(
    mat data,
    colvec theta,
    string family,  // UNUSED
    const colvec order  // UNUSED
) {
  return as<colvec>(cost_gradient(data, theta));
}

CostHessian::CostHessian(Function cost_hessian) :
    cost_hessian(cost_hessian) {}

mat CostHessian::operator()(
    mat data,
    colvec theta,
    string family,  // UNUSED
    double min_prob,  // UNUSED
    const colvec order  // UNUSED
) {
  return as<mat>(cost_hessian(data, theta));
}

}  // namespace fastcpd::wrappers
