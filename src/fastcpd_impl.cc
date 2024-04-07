#include "fastcpd_classes.h"
#include "fastcpd_impl.h"

List fastcpd_impl(
    mat data,
    double beta,
    const string cost_adjustment,
    const unsigned int d,
    const int segment_count,
    const double trim,
    const double momentum_coef,
    Nullable<Function> k,
    const string family,
    const double epsilon,
    const int p,
    const colvec order,
    Nullable<Function> cost,
    Nullable<Function> cost_gradient,
    Nullable<Function> cost_hessian,
    const bool cp_only,
    const double vanilla_percentage,
    const bool warm_start,
    colvec lower,
    colvec upper,
    colvec line_search,
    const mat variance_estimate,
    const unsigned int p_response,
    const double pruning_coef,
    const string r_clock,
    const bool r_progress
) {
  fastcpd::classes::Fastcpd fastcpd_class(
    beta, cost, cost_adjustment, cost_gradient, cost_hessian, cp_only, d, data,
    epsilon, family, k, line_search, lower, momentum_coef, order,
    p, p_response, pruning_coef, r_clock, r_progress, segment_count, trim,
    upper, vanilla_percentage, variance_estimate, warm_start
  );
  return fastcpd_class.run();
}
