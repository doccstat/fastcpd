#include "fastcpd_class.h"
#include "fastcpd_test.h"

using ::fastcpd::classes::Fastcpd;

namespace fastcpd::test {

colvec FastcpdTest::get_gradient_arma(
  const mat& data,
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) {
  Fastcpd fastcpd_class(
    /* beta */ 0,
    /* cost */ R_NilValue,
    /* cost_adjustment */ "MBIC",
    /* cost_gradient */ R_NilValue,
    /* cost_hessian */ R_NilValue,
    /* cp_only */ true,
    /* d */ 0,
    /* data */ data,
    /* epsilon */ 0.0,
    /* family */ "arma",
    /* k */ R_NilValue,
    /* line_search */ colvec(),
    /* lower */ colvec(),
    /* momentum_coef */ 0.0,
    /* order */ colvec({3, 2}),
    /* p */ 0,
    /* p_response */ 0,
    /* pruning_coef */ 0,
    /* r_clock */ "",
    /* r_progress */ false,
    /* segment_count */ 0,
    /* trim */ 0,
    /* upper */ colvec(),
    /* vanilla_percentage */ 0.0,
    /* variance_estimate */ mat(),
    /* warm_start */ false
  );
  return fastcpd_class.get_gradient_arma(segment_start, segment_end, theta);
}

mat FastcpdTest::get_hessian_arma(
  const mat& data,
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) {
  Fastcpd fastcpd_class(
    /* beta */ 0,
    /* cost */ R_NilValue,
    /* cost_adjustment */ "MBIC",
    /* cost_gradient */ R_NilValue,
    /* cost_hessian */ R_NilValue,
    /* cp_only */ true,
    /* d */ 0,
    /* data */ data,
    /* epsilon */ 0.0,
    /* family */ "arma",
    /* k */ R_NilValue,
    /* line_search */ colvec(),
    /* lower */ colvec(),
    /* momentum_coef */ 0.0,
    /* order */ colvec({3, 2}),
    /* p */ 0,
    /* p_response */ 0,
    /* pruning_coef */ 0,
    /* r_clock */ "",
    /* r_progress */ false,
    /* segment_count */ 0,
    /* trim */ 0,
    /* upper */ colvec(),
    /* vanilla_percentage */ 0.0,
    /* variance_estimate */ mat(),
    /* warm_start */ false
  );
  return fastcpd_class.get_hessian_arma(segment_start, segment_end, theta);
}

mat FastcpdTest::get_hessian_binomial(
  const mat& data,
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) {
  Fastcpd fastcpd_class(
    /* beta */ 0,
    /* cost */ R_NilValue,
    /* cost_adjustment */ "MBIC",
    /* cost_gradient */ R_NilValue,
    /* cost_hessian */ R_NilValue,
    /* cp_only */ true,
    /* d */ 0,
    /* data */ data,
    /* epsilon */ 0.0,
    /* family */ "binomial",
    /* k */ R_NilValue,
    /* line_search */ colvec(),
    /* lower */ colvec(),
    /* momentum_coef */ 0.0,
    /* order */ colvec(),
    /* p */ 0,
    /* p_response */ 0,
    /* pruning_coef */ 0,
    /* r_clock */ "",
    /* r_progress */ false,
    /* segment_count */ 0,
    /* trim */ 0,
    /* upper */ colvec(),
    /* vanilla_percentage */ 0.0,
    /* variance_estimate */ mat(),
    /* warm_start */ false
  );
  return fastcpd_class.get_hessian_binomial(segment_start, segment_end, theta);
}

mat FastcpdTest::get_hessian_poisson(
  const mat& data,
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) {
  Fastcpd fastcpd_class(
    /* beta */ 0,
    /* cost */ R_NilValue,
    /* cost_adjustment */ "MBIC",
    /* cost_gradient */ R_NilValue,
    /* cost_hessian */ R_NilValue,
    /* cp_only */ true,
    /* d */ 0,
    /* data */ data,
    /* epsilon */ 0.0,
    /* family */ "poisson",
    /* k */ R_NilValue,
    /* line_search */ colvec(),
    /* lower */ colvec(),
    /* momentum_coef */ 0.0,
    /* order */ colvec(),
    /* p */ 0,
    /* p_response */ 0,
    /* pruning_coef */ 0,
    /* r_clock */ "",
    /* r_progress */ false,
    /* segment_count */ 0,
    /* trim */ 0,
    /* upper */ colvec(),
    /* vanilla_percentage */ 0.0,
    /* variance_estimate */ mat(),
    /* warm_start */ false
  );
  return fastcpd_class.get_hessian_poisson(segment_start, segment_end, theta);
}

double FastcpdTest::get_nll_sen(
    const mat& data,
    const unsigned int segment_start,
    const unsigned int segment_end,
    colvec theta,
    double lambda
) {
  Fastcpd fastcpd_class(
    /* beta */ 0,
    /* cost */ R_NilValue,
    /* cost_adjustment */ "MBIC",
    /* cost_gradient */ R_NilValue,
    /* cost_hessian */ R_NilValue,
    /* cp_only */ true,
    /* d */ 0,
    /* data */ data,
    /* epsilon */ 0.0,
    /* family */ "arma",
    /* k */ R_NilValue,
    /* line_search */ colvec(),
    /* lower */ colvec(),
    /* momentum_coef */ 0.0,
    /* order */ colvec({3, 2}),
    /* p */ 0,
    /* p_response */ 0,
    /* pruning_coef */ 0,
    /* r_clock */ "",
    /* r_progress */ false,
    /* segment_count */ 0,
    /* trim */ 0,
    /* upper */ colvec(),
    /* vanilla_percentage */ 0.0,
    /* variance_estimate */ mat(),
    /* warm_start */ false
  );
  return (fastcpd_class.*fastcpd_class.get_nll_sen)(
    segment_start, segment_end, theta, lambda
  );
}

CostResult FastcpdTest::get_nll_pelt(
  const mat& data,
  const unsigned int segment_start,
  const unsigned int segment_end,
  const double lambda,
  const bool cv,
  const Nullable<colvec>& start
) {
  Fastcpd fastcpd_class(
    /* beta */ 0,
    /* cost */ R_NilValue,
    /* cost_adjustment */ "MBIC",
    /* cost_gradient */ R_NilValue,
    /* cost_hessian */ R_NilValue,
    /* cp_only */ true,
    /* d */ 0,
    /* data */ data,
    /* epsilon */ 0.0,
    /* family */ "arma",
    /* k */ R_NilValue,
    /* line_search */ colvec(),
    /* lower */ colvec(),
    /* momentum_coef */ 0.0,
    /* order */ colvec({3, 2}),
    /* p */ 0,
    /* p_response */ 0,
    /* pruning_coef */ 0,
    /* r_clock */ "",
    /* r_progress */ false,
    /* segment_count */ 0,
    /* trim */ 0,
    /* upper */ colvec(),
    /* vanilla_percentage */ 0.0,
    /* variance_estimate */ mat(),
    /* warm_start */ false
  );
  return (fastcpd_class.*fastcpd_class.get_nll_pelt)(
    segment_start, segment_end, lambda, cv, start
  );
}

mat FastcpdTest::update_theta_sum(
  const unsigned int col,
  colvec old_theta_sum,
  colvec new_theta_sum
) {
  Fastcpd fastcpd_class(
    /* beta */ 0,
    /* cost */ R_NilValue,
    /* cost_adjustment */ "MBIC",
    /* cost_gradient */ R_NilValue,
    /* cost_hessian */ R_NilValue,
    /* cp_only */ true,
    /* d */ 0,
    /* data */ mat(),
    /* epsilon */ 0.0,
    /* family */ "gaussian",
    /* k */ R_NilValue,
    /* line_search */ colvec(),
    /* lower */ colvec(),
    /* momentum_coef */ 0.0,
    /* order */ colvec(),
    /* p */ 3,
    /* p_response */ 0,
    /* pruning_coef */ 0,
    /* r_clock */ "",
    /* r_progress */ false,
    /* segment_count */ 0,
    /* trim */ 0,
    /* upper */ colvec(),
    /* vanilla_percentage */ 0.0,
    /* variance_estimate */ mat(),
    /* warm_start */ false
  );
  fastcpd_class.create_theta_sum(0, old_theta_sum);
  fastcpd_class.update_theta_sum(0, new_theta_sum);
  return fastcpd_class.theta_sum;
}

}  // namespace fastcpd::test
