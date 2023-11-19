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
  return theta.isNull()? cost(data) : cost(data, theta);  // # nocov
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
