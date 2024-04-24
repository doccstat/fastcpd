#ifndef FASTCPD_WRAPPER_H_
#define FASTCPD_WRAPPER_H_

#include "fastcpd_types.h"

namespace fastcpd::classes {

struct ColMat {
  mat data;

  operator colvec() const;
  operator mat() const;
  operator rowvec() const;
};

struct CostResult {
  ColMat par;
  ColMat residuals;
  double value;
};

struct CostFunction {
  const Function cost;
  const mat& data;

  CostFunction(const Function& cost, const mat& data);
  CostResult operator() (
      const unsigned int segment_start,
      const unsigned int segment_end,
      const Nullable<colvec>& theta,
      const double lambda,  // UNUSED
      const bool cv,  // UNUSED
      const Nullable<colvec>& start  // UNUSED
  ) const;
};

struct CostGradient {
  const Function cost_gradient;
  const mat& data;

  CostGradient(const Function& cost_gradient, const mat& data);
  const colvec operator() (
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  ) const;
};

struct CostHessian {
  const Function cost_hessian;
  const mat& data;

  CostHessian(const Function& cost_hessian, const mat& data);
  const mat operator() (
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  ) const;
};

}  // namespace fastcpd::classes

#endif  // FASTCPD_WRAPPER_H_
