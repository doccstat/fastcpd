#ifndef FASTCPD_FAMILIES_MEANVARIANCE_H_
#define FASTCPD_FAMILIES_MEANVARIANCE_H_

#include "fastcpd_family.h"

namespace fastcpd::families {

struct MeanvarianceFamily : BaseFamily {
  static constexpr const char* name = "meanvariance";
  static constexpr bool has_optimized_run = false;
  static constexpr bool is_pelt_only = true;

  // Data layout: TRANSPOSED — shape (p+p²) × (n+1) instead of (n+1) × (p+p²).
  // Same reasoning as VarianceFamily: .row(t) in column-major strides by n+1
  // doubles between each element; .col(t) after transpose is contiguous.
  // data_c_n_rows_ == p+p² after the transpose.
  static arma::mat CreateDataC(arma::mat const& data,
                               arma::mat const& variance_estimate,
                               unsigned int const p_response = 0) {
    arma::mat data_crossprod(data.n_cols * data.n_cols, data.n_rows);
    for (unsigned int i = 0; i < data.n_rows; i++) {
      data_crossprod.col(i) = arma::vectorise(data.row(i).t() * data.row(i));
    }
    arma::mat data_c = arma::cumsum(arma::join_rows(data, data_crossprod.t()));
    data_c = arma::join_cols(arma::zeros<arma::rowvec>(data_c.n_cols), data_c);
    return data_c.t();  // (p+p²) × (n+1): time is the column index
  }

  static unsigned int GetDataNDims(arma::mat const& data) {
    return data.n_cols;
  }

  template <typename Solver>
  static void GetNllPelt(Solver* solver, unsigned int const segment_start,
                         unsigned int const segment_end, bool const cv,
                         std::optional<arma::colvec> const& start) {
    arma::mat data_segment = solver->data_.rows(segment_start, segment_end);
    arma::mat cov_matrix = arma::cov(data_segment);
    solver->result_coefficients_ =
        arma::zeros<arma::colvec>(solver->parameters_count_);
    solver->result_coefficients_.rows(0, solver->data_n_dims_ - 1) =
        arma::mean(data_segment, 0).t();
    solver->result_coefficients_.rows(solver->data_n_dims_,
                                      solver->parameters_count_ - 1) =
        cov_matrix.as_col();
    solver->result_residuals_ =
        data_segment.each_row() - arma::mean(data_segment, 0);
    solver->result_residuals_.each_row() /= arma::sqrt(cov_matrix.diag()).t();
    GetNllPeltValue(solver, segment_start, segment_end, cv, start);
  }

  template <typename Solver>
  static void GetNllPeltValue(Solver* solver, unsigned int const segment_start,
                              unsigned int const segment_end, bool const cv,
                              std::optional<arma::colvec> const& start) {
    solver->result_value_ = GetNllPeltValueFast<-1>(solver, segment_start, segment_end);
  }

  // With the transposed layout data_c_.col(t) returns (p+p²) contiguous
  // doubles.  Indices 0..p-1 hold Σx; indices p + j2*p + j1 hold Σ(x[j1]*x[j2])
  // (column-major vectorise order of the outer product x*x.t()).
  //
  // When kNDims == 1 (p == 1), data_c_ is 2×(n+1): [Σx, Σx²] per col.
  // The det of a 1×1 sample-covariance matrix is just the scalar variance.
  //
  // General path (p > 1): build the p×p sample-covariance matrix directly
  // from raw pointers — no arma::colvec temporaries, no heap alloc per
  // candidate.  Previously three heap objects were created per call (data_diff,
  // mean_diff, mean_diff*mean_diff.t()), which showed up as 2+ s system time
  // for n=1e6 examples.
  template <int kNDims, typename Solver>
  static double GetNllPeltValueFast(Solver* solver, unsigned int const segment_start,
                                    unsigned int const segment_end) {
    unsigned int const segment_length = segment_end - segment_start + 1;
    double const n = static_cast<double>(segment_length);
    if constexpr (kNDims == 1) {
      // data_c_ is (2) × (n+1), stride = 2: [cumsum_x, cumsum_x²] per col.
      double const* const end_ptr =
          solver->data_c_ptr_ + static_cast<std::size_t>(segment_end + 1) * 2;
      double const* const start_ptr =
          solver->data_c_ptr_ + static_cast<std::size_t>(segment_start) * 2;
      double const sum_x = end_ptr[0] - start_ptr[0];
      double const sum_xx = end_ptr[1] - start_ptr[1];
      double const det_value =
          (sum_xx - sum_x * sum_x / n) / n;
      return std::log(det_value > 0.0 ? det_value : 1e-10) * n / 2.0;
    }
    // General p > 1: raw pointer path.
    // ep[j]         = cumsum_x[j]           (j = 0..p-1)
    // ep[p+j2*p+j1] = cumsum of x[j1]*x[j2] (column-major vectorise order)
    unsigned int const p = solver->data_n_dims_;
    double const* const ep =
        solver->data_c_ptr_ +
        static_cast<std::size_t>(segment_end + 1) * solver->data_c_n_rows_;
    double const* const sp =
        solver->data_c_ptr_ +
        static_cast<std::size_t>(segment_start) * solver->data_c_n_rows_;
    arma::mat cov_mat(p, p);
    for (unsigned int j1 = 0; j1 < p; j1++) {
      double const sum_j1 = ep[j1] - sp[j1];
      for (unsigned int j2 = 0; j2 < p; j2++) {
        double const sum_j2 = ep[j2] - sp[j2];
        cov_mat(j1, j2) =
            (ep[p + j2 * p + j1] - sp[p + j2 * p + j1] -
             sum_j1 * sum_j2 / n) / n;
      }
    }
    double det_value = arma::det(cov_mat);
    if (det_value <= 0) det_value = 1e-10;
    return std::log(det_value) * n / 2.0;
  }

  // Prefetch all p+p² values for candidate s — contiguous after the transpose.
  // Loop over 64-byte cache lines: no branch, handles any p.
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

#endif  // FASTCPD_FAMILIES_MEANVARIANCE_H_
