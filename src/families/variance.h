#ifndef FASTCPD_FAMILIES_VARIANCE_H_
#define FASTCPD_FAMILIES_VARIANCE_H_

#include "fastcpd_family.h"

#include <vector>

namespace fastcpd::families {

struct VarianceFamily : BaseFamily {
  static constexpr const char* name = "variance";
  static constexpr bool has_optimized_run = false;
  static constexpr bool is_pelt_only = true;

  // Data layout: TRANSPOSED — shape (p²) × (n+1) instead of (n+1) × p².
  // data_c_.row(t) in column-major strides by n+1 doubles between each of the
  // p² elements (≈ 8 GB apart for n = 10^9).  After transposing, data_c_.col(t)
  // returns p² contiguous doubles — one cache line for p ≤ 2, eight for p = 8.
  // data_c_n_rows_ == p² after the transpose and serves as the prefetch stride.
  static arma::mat CreateDataC(arma::mat const& data,
                               arma::mat const& variance_estimate,
                               unsigned int const p_response = 0) {
    arma::uword const n = data.n_rows;
    arma::uword const p = data.n_cols;
    arma::uword const total = p * p;
    arma::rowvec const means = arma::mean(data, 0);
    arma::mat out(total, n + 1, arma::fill::none);
    out.col(0).zeros();

    if (p == 1) {
      double const mean = means[0];
      double prefix = 0.0;
      double const* const x = data.colptr(0);
      for (arma::uword i = 0; i < n; i++) {
        double const centered = x[i] - mean;
        prefix += centered * centered;
        out(0, i + 1) = prefix;
      }
      return out;  // 1 × (n+1): time is the column index
    }

    std::vector<double> centered(p);
    std::vector<double> prefix(total, 0.0);
    for (arma::uword i = 0; i < n; i++) {
      for (arma::uword j = 0; j < p; j++) {
        centered[j] = data(i, j) - means[j];
      }
      arma::uword idx = 0;
      for (arma::uword j2 = 0; j2 < p; j2++) {
        double const x_j2 = centered[j2];
        for (arma::uword j1 = 0; j1 < p; j1++) {
          prefix[idx] += centered[j1] * x_j2;
          out(idx, i + 1) = prefix[idx];
          idx++;
        }
      }
    }
    return out;  // (p²) × (n+1): time is the column index
  }

  static unsigned int GetDataNDims(arma::mat const& data) {
    return data.n_cols;
  }

  template <typename Solver>
  static void GetNllPelt(Solver* solver, unsigned int const segment_start,
                         unsigned int const segment_end, bool const cv,
                         std::optional<arma::colvec> const& start) {
    arma::mat data_segment = solver->data_.rows(segment_start, segment_end);
    arma::mat covar = arma::cov(data_segment);
    solver->result_coefficients_ = covar.as_col();
    solver->result_residuals_ =
        data_segment.each_row() / arma::sqrt(covar.diag()).t();
    GetNllPeltValue(solver, segment_start, segment_end, cv, start);
  }

  template <typename Solver>
  static void GetNllPeltValue(Solver* solver, unsigned int const segment_start,
                              unsigned int const segment_end, bool const cv,
                              std::optional<arma::colvec> const& start) {
    solver->result_value_ = GetNllPeltValueFast<-1>(solver, segment_start, segment_end);
  }

  // With the transposed layout data_c_.col(t) returns p² contiguous doubles.
  // Indices are in column-major vectorise order: j2*p+j1 holds cumsum of
  // (x-mean)[j1] * (x-mean)[j2]  (data is pre-centred in CreateDataC).
  //
  // When kNDims == 1 (p == 1), data_c_ is 1×(n+1): direct scalar path.
  //
  // General path (p > 1): build the p×p covariance matrix directly from raw
  // pointers — no arma::colvec subtraction, no arma::reshape, no / scalar
  // temporary — eliminating 2–3 heap allocs per candidate vs the old
  // arma::det(arma::reshape(col(...) - col(...), p, p) / n) expression.
  template <int kNDims, typename Solver>
  static double GetNllPeltValueFast(Solver* solver, unsigned int const segment_start,
                                    unsigned int const segment_end) {
    if constexpr (kNDims == 1) {
      // data_c_ is (1) × (n+1), stride = 1.
      unsigned int const segment_length = segment_end - segment_start + 1;
      double const sum_sq =
          solver->data_c_ptr_[segment_end + 1] - solver->data_c_ptr_[segment_start];
      double const seg_var = sum_sq / static_cast<double>(segment_length);
      return std::log(seg_var > 0.0 ? seg_var : 1e-10) *
             static_cast<double>(segment_length) / 2.0;
    }
    // Adjust bounds for short segments (fewer than p observations).
    unsigned int approx_start = segment_start, approx_end = segment_end;
    if (approx_end - approx_start + 1 < solver->data_n_dims_) {
      if (segment_end < solver->data_n_rows_ - solver->data_n_dims_) {
        approx_end = segment_end + solver->data_n_dims_;
      } else {
        approx_end = solver->data_n_rows_ - 1;
      }
      approx_start = approx_end - solver->data_n_dims_;
    }
    unsigned int const p = solver->data_n_dims_;
    double const n = static_cast<double>(approx_end - approx_start + 1);
    double const* const ep =
        solver->data_c_ptr_ +
        static_cast<std::size_t>(approx_end + 1) * solver->data_c_n_rows_;
    double const* const sp =
        solver->data_c_ptr_ +
        static_cast<std::size_t>(approx_start) * solver->data_c_n_rows_;
    arma::mat cov_mat(p, p);
    for (unsigned int j2 = 0; j2 < p; j2++) {
      for (unsigned int j1 = 0; j1 < p; j1++) {
        cov_mat(j1, j2) = (ep[j2 * p + j1] - sp[j2 * p + j1]) / n;
      }
    }
    double const det_value = arma::det(cov_mat);
    return std::log(det_value) * n / 2.0;
  }

  // Prefetch all p² values for candidate s — they are contiguous after the
  // transpose.  Loop over 64-byte cache lines: no branch, handles any p.
  template <typename Solver>
  static void PrefetchCandidate(Solver* solver, unsigned int const s) {
    double const* const ptr =
        solver->data_c_ptr_ + static_cast<std::size_t>(s) * solver->data_c_n_rows_;
    for (unsigned int b = 0; b < solver->data_c_n_rows_; b += 8) {
      absl::PrefetchToLocalCache(ptr + b);
    }
  }
};

}  // namespace fastcpd::families

#endif  // FASTCPD_FAMILIES_VARIANCE_H_
