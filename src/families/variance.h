#ifndef FASTCPD_FAMILIES_VARIANCE_H_
#define FASTCPD_FAMILIES_VARIANCE_H_

#include "fastcpd_family.h"

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
    arma::mat data_c = data;
    data_c.each_row() -= arma::mean(data_c, 0);
    arma::mat data_crossprod(data.n_cols * data.n_cols, data.n_rows);
    for (unsigned int i = 0; i < data.n_rows; i++) {
      data_crossprod.col(i) =
          arma::vectorise(data_c.row(i).t() * data_c.row(i));
    }
    data_c = arma::cumsum(data_crossprod.t());
    data_c = arma::join_cols(arma::zeros<arma::rowvec>(data_c.n_cols), data_c);
    return data_c.t();  // (p²) × (n+1): time is the column index
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
#if defined(__GNUC__) || defined(__clang__)
    double const* const ptr =
        solver->data_c_ptr_ + static_cast<std::size_t>(s) * solver->data_c_n_rows_;
    for (unsigned int b = 0; b < solver->data_c_n_rows_; b += 8) {
      __builtin_prefetch(ptr + b, 0, 1);
    }
#endif
  }
};

}  // namespace fastcpd::families

#endif  // FASTCPD_FAMILIES_VARIANCE_H_
