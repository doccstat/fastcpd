#ifndef FASTCPD_TEST_H_
#define FASTCPD_TEST_H_

#include "fastcpd_wrapper.h"

using ::fastcpd::classes::CostResult;
using ::Rcpp::Nullable;

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

  static double get_nll_sen(
    const mat& data,
    const unsigned int segment_start,
    const unsigned int segment_end,
    colvec theta,
    double lambda
  );

  static CostResult get_nll_pelt(
    const mat& data,
    const unsigned int segment_start,
    const unsigned int segment_end,
    const double lambda,
    const bool cv,
    const Nullable<colvec>& start
  );

  static mat update_theta_sum(
    const unsigned int col,
    colvec old_theta_sum,
    colvec new_theta_sum
  );
};

}  // namespace fastcpd::test

#endif  // FASTCPD_TEST_H_
