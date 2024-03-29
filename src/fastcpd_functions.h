#ifndef FASTCPD_FUNCTIONS_H_
#define FASTCPD_FUNCTIONS_H_

#include "fastcpd_types.h"
#include "fastcpd_classes.h"

using ::fastcpd::classes::CostResult;

namespace fastcpd::functions {

CostResult negative_log_likelihood_arma(
  const mat data,
  const colvec order
);

CostResult negative_log_likelihood_glm(
  const mat data,
  Nullable<colvec> start,
  const std::string family
);

CostResult negative_log_likelihood_lasso_cv(const mat data);

CostResult negative_log_likelihood_lasso_wo_cv(
  const mat data,
  const double lambda
);

CostResult negative_log_likelihood_mean(
  const mat data,
  const mat variance_estimate
);

CostResult negative_log_likelihood_meanvariance(
  const mat data,
  const double epsilon
);

CostResult negative_log_likelihood_mgaussian(
  const mat data,
  const unsigned int p_response,
  const mat variance_estimate
);

CostResult negative_log_likelihood_variance(
  const mat data,
  const rowvec variance_data_mean
);

}  // namespace fastcpd::functions

#endif  // FASTCPD_FUNCTIONS_H_
