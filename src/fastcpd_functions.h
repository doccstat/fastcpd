#ifndef FASTCPD_FUNCTIONS_H_
#define FASTCPD_FUNCTIONS_H_

#include "fastcpd_types.h"
#include "fastcpd_classes.h"

using ::fastcpd::classes::CostResult;

namespace fastcpd::functions {

CostResult negative_log_likelihood_lasso_cv(const mat data);

CostResult negative_log_likelihood_mean(
  const mat data,
  const mat variance_estimate
);

CostResult negative_log_likelihood_meanvariance(
  const mat data,
  const double epsilon
);

}  // namespace fastcpd::functions

#endif  // FASTCPD_FUNCTIONS_H_
