#ifndef FASTCPD_FAMILIES_EXPONENTIAL_H_
#define FASTCPD_FAMILIES_EXPONENTIAL_H_

#include "fastcpd_family.h"

namespace fastcpd::families {

// Mean-change detection under exponentially distributed noise: each segment
// is modeled as i.i.d. Exponential(rate = lambda), and the unknown rate
// changes across segments (cf. `changepoint::cpt.meanvar` with
// `test.stat = "Exponential"`).
//
// For a segment of length n with observations x_1..x_n, the rate MLE is
//   lambda_hat = n / sum(x) = 1 / mean(x)
// and the resulting per-segment negative log-likelihood has the closed form
//   NLL = n * (log(mean(x)) + 1)
// which — like `mean` and `variance` — supports an O(1) cumulative-sum
// lookup per PELT candidate.
struct ExponentialFamily : BaseFamily {
  static constexpr const char* name = "exponential";
  static constexpr bool has_optimized_run = false;
  static constexpr bool is_pelt_only = true;

  // This family is intrinsically univariate (a single rate parameter), so
  // `data_c_` is just the cumulative sum of the raw observations, stored
  // transposed (1 x (n+1)) to match the cache-friendly layout convention used
  // by `mean`/`variance`/`meanvariance`.
  static arma::mat CreateDataC(arma::mat const& data,
                               arma::mat const& variance_estimate,
                               unsigned int const p_response = 0) {
    arma::mat data_c = arma::cumsum(data);
    data_c = arma::join_cols(arma::zeros<arma::rowvec>(data_c.n_cols), data_c);
    return data_c.t();  // 1 x (n+1)
  }

  static unsigned int GetDataNDims(arma::mat const& data) {
    return data.n_cols;
  }

  template <typename Solver>
  static void GetNllPelt(Solver* solver, unsigned int const segment_start,
                         unsigned int const segment_end, bool const cv,
                         std::optional<arma::colvec> const& start) {
    arma::mat const data_segment =
        solver->data_.rows(segment_start, segment_end);
    arma::rowvec const segment_mean = arma::mean(data_segment, 0);
    solver->result_coefficients_ = 1.0 / segment_mean.t();
    solver->result_residuals_ =
        data_segment.each_row() / segment_mean -
        arma::ones<arma::mat>(data_segment.n_rows, data_segment.n_cols);
    GetNllPeltValue(solver, segment_start, segment_end, cv, start);
  }

  template <typename Solver>
  static void GetNllPeltValue(Solver* solver, unsigned int const segment_start,
                              unsigned int const segment_end, bool const cv,
                              std::optional<arma::colvec> const& start) {
    solver->result_value_ =
        GetNllPeltValueFast<-1>(solver, segment_start, segment_end);
  }

  // O(1) closed-form NLL via cumulative-sum lookup:
  //   segment_sum  = cumsum[end+1] - cumsum[start]
  //   segment_mean = segment_sum / segment_length
  //   NLL          = segment_length * (log(segment_mean) + 1)
  //
  // The family has exactly one parameter (the rate), so unlike `mean` there
  // is no general-`p` branch to specialize away from: `kNDims` is accepted
  // only to satisfy `GetCostValuePelt`'s dispatch signature and is otherwise
  // ignored — the scalar formula below is always the right one.
  template <int kNDims, typename Solver>
  static double GetNllPeltValueFast(Solver* solver,
                                    unsigned int const segment_start,
                                    unsigned int const segment_end) {
    double const* const data_c = solver->data_c_ptr_;
    unsigned int const stride = solver->data_c_n_rows_;  // == 1
    double const segment_sum =
        data_c[static_cast<std::size_t>(segment_end + 1) * stride] -
        data_c[static_cast<std::size_t>(segment_start) * stride];
    unsigned int const segment_length = segment_end - segment_start + 1;
    double const segment_mean = segment_sum / segment_length;
    return segment_length * (std::log(segment_mean) + 1.0);
  }

  template <typename Solver>
  static void PrefetchCandidate(Solver* solver, unsigned int const s) {
    absl::PrefetchToLocalCache(solver->data_c_ptr_ + s);
  }
};

}  // namespace fastcpd::families

#endif  // FASTCPD_FAMILIES_EXPONENTIAL_H_
