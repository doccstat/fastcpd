#ifndef FASTCPD_FAMILY_H_
#define FASTCPD_FAMILY_H_

#include <armadillo>
#include <cstring>
#include <optional>
#include <string>

namespace fastcpd {
namespace classes {

enum class CostAdjustment { kBIC, kMBIC, kMDL };

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
class Fastcpd;
}

namespace families {

struct BaseFamily {
  // True for families that never use sequential gradient descent (SEN).
  // Eliminates all SEN-related code paths at compile time.
  static constexpr bool is_pelt_only = false;
  // True for families that support warm-start initialisation (binomial/poisson).
  static constexpr bool supports_warm_start = false;
  // True only for CustomFamily (user-supplied cost functions).
  static constexpr bool is_custom = false;
  // True for families that track per-segment error standard deviation
  // (gaussian / lasso), used for MBIC penalty and lasso beta scaling.
  static constexpr bool has_error_std_dev = false;
  // True only for LassoFamily (needs soft-thresholding in the SEN step).
  static constexpr bool is_lasso = false;

  static arma::mat CreateDataC(arma::mat const& data,
                               arma::mat const& variance_estimate,
                               unsigned int const p_response = 0) {
    return data;
  }

  static unsigned int GetDataNDims(arma::mat const& data) {
    return data.n_cols;
  }

  // Default: initialize SEN parameters from first segment statistics.
  template <typename Solver>
  static void CreateSenParameters(Solver* solver) {
    solver->coefficients_.col(0) = solver->segment_coefficients_.row(0).t();
    solver->coefficients_sum_.col(0) = solver->segment_coefficients_.row(0).t();
    solver->hessian_.slice(0) =
        solver->epsilon_in_hessian_ *
        arma::eye<arma::mat>(solver->parameters_count_,
                             solver->parameters_count_);
  }

  // Default: seed the newest candidate slot with the current segment's
  // coefficients and an identity-scaled Hessian.
  template <typename Solver>
  static void UpdateSenParameters(Solver* solver) {
    int const segment_index =
        arma::index_max(arma::find(solver->segment_indices_ <= solver->t - 1));
    arma::colvec coef_add = solver->segment_coefficients_.row(segment_index).t();
    arma::mat hessian_new =
        solver->epsilon_in_hessian_ *
        arma::eye<arma::mat>(solver->parameters_count_,
                             solver->parameters_count_);
    std::memcpy(solver->coefficients_.colptr(solver->pruned_set_size_ - 1),
                coef_add.memptr(),
                sizeof(double) * solver->parameters_count_);
    std::memcpy(solver->coefficients_sum_.colptr(solver->pruned_set_size_ - 1),
                coef_add.memptr(),
                sizeof(double) * solver->parameters_count_);
    std::memcpy(
        solver->hessian_.slice(solver->pruned_set_size_ - 1).memptr(),
        hessian_new.memptr(),
        sizeof(double) * solver->parameters_count_ * solver->parameters_count_);
  }

  // Default: PELT-only families that do not implement SEN return 0.
  // This is only reached when SEN is inadvertently triggered (e.g. during
  // initialisation with vanilla_percentage != 1); the optimised Run() path
  // overwrites any result produced by this call.
  template <typename Solver>
  static double GetNllSen(Solver* solver, unsigned int const segment_start,
                          unsigned int const segment_end,
                          arma::colvec const& theta) {
    return 0.0;
  }

  template <typename Solver>
  static void GetNllPeltValue(Solver* solver,
                              unsigned int const segment_start,
                              unsigned int const segment_end, bool const cv,
                              std::optional<arma::colvec> const& start) {
    // Default implementation
  }

  // Called from UpdateStep before evaluating candidate s, with a lookahead
  // of kPrefetchDist steps, to hide DRAM latency for data_c_ accesses.
  // Default is a no-op; PELT families with direct data_c_ptr_ access override.
  template <typename Solver>
  static void PrefetchCandidate(Solver* solver, unsigned int const s) {}

  template <typename Solver>
  static arma::colvec GetGradient(Solver* solver,
                                  unsigned int const segment_start,
                                  unsigned int const segment_end,
                                  arma::colvec const& theta) {
    return arma::ones(theta.n_elem);
  }

  template <typename Solver>
  static arma::mat GetHessian(Solver* solver,
                              unsigned int const segment_start,
                              unsigned int const segment_end,
                              arma::colvec const& theta) {
    return arma::eye(theta.n_elem, theta.n_elem);
  }
};

}  // namespace families
}  // namespace fastcpd

#endif  // FASTCPD_FAMILY_H_
