#ifndef FASTCPD_TEST_H_
#define FASTCPD_TEST_H_

#include "fastcpd_types.h"
#include "fastcpd_wrapper.h"

using ::fastcpd::classes::CostResult;

namespace fastcpd::test {

class FastcpdTest {
 public:
  static colvec get_gradient_arma(
    const mat& data,
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  static mat get_hessian_arma(
    const mat& data,
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  static mat get_hessian_binomial(
    const mat& data,
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  static mat get_hessian_poisson(
    const mat& data,
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  static double get_nll_wo_cv(
    const mat& data,
    const unsigned int segment_start,
    const unsigned int segment_end,
    colvec theta,
    double lambda
  );

  static CostResult get_nll_wo_theta(
    const mat& data,
    const unsigned int segment_start,
    const unsigned int segment_end,
    double lambda,
    bool cv,
    Nullable<colvec> start
  );

  static mat update_theta_sum(
    const unsigned int col,
    colvec old_theta_sum,
    colvec new_theta_sum
  );
};

}  // namespace fastcpd::test

#endif  // FASTCPD_TEST_H_
