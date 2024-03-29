#ifndef FASTCPD_FUNCTIONS_H_
#define FASTCPD_FUNCTIONS_H_

#include "fastcpd_types.h"
#include "fastcpd_classes.h"

namespace fastcpd::functions {

fastcpd::classes::CostResult negative_log_likelihood_mean(
  const mat data,
  const mat variance_estimate
);

}  // namespace fastcpd::functions

#endif  // FASTCPD_FUNCTIONS_H_
