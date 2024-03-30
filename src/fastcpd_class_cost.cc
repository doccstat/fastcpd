#include "fastcpd_classes.h"

namespace fastcpd::classes {

ColMat::operator colvec() const {
  // TODO(doccstat): Add a warning if the matrix has more than one column.
  return data.as_col();
}

ColMat::operator mat() const {
  return data;
}

ColMat::operator rowvec() const {
  // TODO(doccstat): Add a warning if the matrix has more than one column.
  return data.as_col().t();
}

CostFunction::CostFunction(Function cost) : cost(cost) {}

CostResult CostFunction::operator()(  // # nocov
    mat data,
    Nullable<colvec> theta,
    double lambda,  // UNUSED
    bool cv,  // UNUSED
    Nullable<colvec> start  // UNUSED
) {
  DEBUG_RCOUT(data.n_rows);
  SEXP value =
    theta.isNull() ? cost(data) : cost(data, as<colvec>(theta));  // # nocov
  return {{colvec()}, {colvec()}, as<double>(value)};  // # nocov
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
