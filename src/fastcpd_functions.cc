#include "fastcpd_functions.h"

using ::fastcpd::classes::CostResult;

namespace fastcpd::functions {

CostResult negative_log_likelihood_mean(const mat data, const mat variance_estimate) {
  rowvec par = mean(data, 0);
  return CostResult{
    par,
    data.each_row() - par,
    data.n_rows / 2.0 * (
      std::log(2.0 * M_PI) * data.n_cols + log_det_sympd(variance_estimate) +
        trace(solve(
          variance_estimate,
          (data.each_row() - par).t() * (data.each_row() - par)
        )) / data.n_rows
    )
  };
}

}  // namespace fastcpd::functions
