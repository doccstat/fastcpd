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

CostFunction::CostFunction(
  const Function& cost,
  const mat& data
) : cost(cost),
    data(data) {}

CostResult CostFunction::operator() (  // # nocov
    const unsigned int segment_start,
    const unsigned int segment_end,
    const Nullable<colvec>& theta,
    const double lambda,  // UNUSED
    const bool cv,  // UNUSED
    const Nullable<colvec>& start  // UNUSED
) const {
  SEXP value = theta.isNull() ?
    cost(data.rows(segment_start, segment_end)) :
    cost(data.rows(segment_start, segment_end), as<colvec>(theta));  // # nocov
  return {{colvec()}, {colvec()}, as<double>(value)};  // # nocov
}

CostGradient::CostGradient(
  const Function& cost_gradient,
  const mat& data
) : cost_gradient(cost_gradient),
    data(data) {}

const colvec CostGradient::operator() (
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) const {
  return as<colvec>(
    cost_gradient(data.rows(segment_start, segment_end), theta)
  );
}

CostHessian::CostHessian(
  const Function& cost_hessian,
  const mat& data
) : cost_hessian(cost_hessian),
    data(data) {}

const mat CostHessian::operator() (
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) const {
  return as<mat>(cost_hessian(data.rows(segment_start, segment_end), theta));
}

}  // namespace fastcpd::classes
