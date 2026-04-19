#ifndef FASTCPD_FAMILIES_MGAUSSIAN_H_
#define FASTCPD_FAMILIES_MGAUSSIAN_H_

#include <cmath>

#include "fastcpd_family.h"

namespace fastcpd::families {

// Multivariate Gaussian (mgaussian) family for VAR models.
// Data layout: first regression_response_count_ columns are responses (y),
// remaining columns are predictors (x).
struct MgaussianFamily : BaseFamily {
  static constexpr const char* name = "mgaussian";
  static constexpr bool has_optimized_run = false;
  static constexpr bool is_pelt_only = true;

  static unsigned int GetDataNDims(arma::mat const& data) {
    return data.n_cols;
  }

  // Data layout: TRANSPOSED — shape (q²+p·q+p²) × (n+1).
  // Indices per time column (column-major vectorise order):
  //   [0 .. q²-1]          YY: Σ y[j1]·y[j2],  j2*q+j1
  //   [q² .. q²+p·q-1]     XY: Σ x[j1]·y[j2],  j2*p+j1
  //   [q²+p·q .. total-1]  XX: Σ x[j1]·x[j2],  j2*p+j1
  // where y = first q cols, x = remaining p cols.
  static arma::mat CreateDataC(arma::mat const& data,
                               arma::mat const& variance_estimate,
                               unsigned int const p_response) {
    unsigned int const n = data.n_rows;
    unsigned int const q = p_response;
    unsigned int const p = data.n_cols - q;  // predictor count
    unsigned int const yy_size = q * q;
    unsigned int const xy_size = p * q;
    unsigned int const xx_size = p * p;
    unsigned int const total = yy_size + xy_size + xx_size;

    // Build the total×n cross-product matrix, then cumsum + transpose.
    // Layout per column i: [YY flat | XY flat | XX flat] in column-major order.
    arma::mat data_crossprod(total, n);
    for (unsigned int i = 0; i < n; i++) {
      unsigned int idx = 0;
      // YY: q×q outer product y*y', column-major: (j1,j2) → j2*q+j1
      for (unsigned int j2 = 0; j2 < q; j2++)
        for (unsigned int j1 = 0; j1 < q; j1++)
          data_crossprod(idx++, i) = data(i, j1) * data(i, j2);
      // XY: p×q outer product x*y', column-major: (j1,j2) → j2*p+j1
      for (unsigned int j2 = 0; j2 < q; j2++)
        for (unsigned int j1 = 0; j1 < p; j1++)
          data_crossprod(idx++, i) = data(i, q + j1) * data(i, j2);
      // XX: p×p outer product x*x', column-major: (j1,j2) → j2*p+j1
      for (unsigned int j2 = 0; j2 < p; j2++)
        for (unsigned int j1 = 0; j1 < p; j1++)
          data_crossprod(idx++, i) = data(i, q + j1) * data(i, q + j2);
    }
    // cumsum along rows (dim=0) of the n×total transpose, then prepend zeros.
    arma::mat cum = arma::cumsum(data_crossprod.t());  // n × total
    cum = arma::join_cols(arma::zeros<arma::rowvec>(total), cum);  // (n+1) × total
    return cum.t();  // total × (n+1): time is the column index
  }

  template <typename Solver>
  static void GetNllPelt(Solver* solver, unsigned int const segment_start,
                         unsigned int const segment_end, bool const cv,
                         std::optional<arma::colvec> const& start) {
    arma::mat const data_segment =
        solver->data_.rows(segment_start, segment_end);
    unsigned int const q = solver->regression_response_count_;
    arma::mat x = data_segment.cols(q, data_segment.n_cols - 1);
    arma::mat y = data_segment.cols(0, q - 1);
    arma::mat x_t_x;
    if (data_segment.n_rows <= data_segment.n_cols - q + 1) {
      x_t_x = arma::eye<arma::mat>(data_segment.n_cols - q,
                                   data_segment.n_cols - q);
    } else {
      x_t_x = x.t() * x;
    }
    arma::mat coefficients = arma::solve(x_t_x, x.t()) * y;
    solver->result_coefficients_ = coefficients.as_col();
    solver->result_residuals_ = y - x * coefficients;
    // Use precomputed log|Σ| to avoid calling log_det_sympd per segment.
    double value =
        static_cast<double>(q) * std::log(2.0 * M_PI) +
        solver->variance_log_det_;
    value *= static_cast<double>(data_segment.n_rows);
    value += arma::trace(arma::solve(solver->variance_estimate_,
                                     solver->result_residuals_.t() *
                                         solver->result_residuals_));
    solver->result_value_ = value / 2.0;
  }

  template <typename Solver>
  static void GetNllPeltValue(Solver* solver, unsigned int const segment_start,
                              unsigned int const segment_end, bool const cv,
                              std::optional<arma::colvec> const& start) {
    solver->result_value_ =
        GetNllPeltValueFast<-1>(solver, segment_start, segment_end);
  }

  // O(p³+q³) per candidate: reads cumulative sums from data_c_,
  // builds XX, XY, YY via pointer subtraction, then solves
  // B = XX⁻¹ XY and computes RSS = YY − XY' B without touching data rows.
  //
  // NLL = n/2·(q·log(2π) + log|Σ|) + trace(Σ⁻¹·RSS)/2
  //
  // Short segment (n ≤ p+1): falls back to GetNllPelt (the Schur complement
  // identity RSS = Y'Y − (X'Y)'(X'X)⁻¹(X'Y) does NOT hold when X'X is
  // replaced by an identity regulariser, so the fast formula would give a
  // different value than GetNllPelt for these degenerate segments).
  template <int kNDims, typename Solver>
  static double GetNllPeltValueFast(Solver* solver,
                                    unsigned int const segment_start,
                                    unsigned int const segment_end) {
    unsigned int const q = solver->regression_response_count_;
    unsigned int const p = solver->data_n_dims_ - q;
    unsigned int const n = segment_end - segment_start + 1;

    // Short-segment fallback: Schur complement requires the actual (X'X)⁻¹.
    if (n <= p + 1) {
      GetNllPelt(solver, segment_start, segment_end, false, std::nullopt);
      return solver->result_value_;
    }

    double const nd = static_cast<double>(n);
    unsigned int const xy_offset = q * q;
    unsigned int const xx_offset = q * q + p * q;

    double const* const ep =
        solver->data_c_ptr_ +
        static_cast<std::size_t>(segment_end + 1) * solver->data_c_n_rows_;
    double const* const sp =
        solver->data_c_ptr_ +
        static_cast<std::size_t>(segment_start) * solver->data_c_n_rows_;

    // Build YY (q×q), XY (p×q), XX (p×p) from cumulative differences.
    arma::mat yy(q, q);
    for (unsigned int j2 = 0; j2 < q; j2++)
      for (unsigned int j1 = 0; j1 < q; j1++)
        yy(j1, j2) = ep[j2 * q + j1] - sp[j2 * q + j1];

    arma::mat xy(p, q);
    for (unsigned int j2 = 0; j2 < q; j2++)
      for (unsigned int j1 = 0; j1 < p; j1++)
        xy(j1, j2) = ep[xy_offset + j2 * p + j1] - sp[xy_offset + j2 * p + j1];

    arma::mat xx(p, p);
    for (unsigned int j2 = 0; j2 < p; j2++)
      for (unsigned int j1 = 0; j1 < p; j1++)
        xx(j1, j2) = ep[xx_offset + j2 * p + j1] - sp[xx_offset + j2 * p + j1];

    // B = (X'X)⁻¹ X'Y  (p×q); RSS = Y'Y − (X'Y)' B  (q×q).
    arma::mat const b = arma::solve(xx, xy);
    arma::mat const rss = yy - xy.t() * b;

    return nd / 2.0 *
               (static_cast<double>(q) * std::log(2.0 * M_PI) +
                solver->variance_log_det_) +
           arma::trace(arma::solve(solver->variance_estimate_, rss)) / 2.0;
  }

  // Prefetch all (q²+p·q+p²) values for candidate s — contiguous after
  // the transpose.  Same cache-line loop as variance/meanvariance.
  template <typename Solver>
  static void PrefetchCandidate(Solver* solver, unsigned int const s) {
#if defined(__GNUC__) || defined(__clang__)
    double const* const ptr =
        solver->data_c_ptr_ +
        static_cast<std::size_t>(s) * solver->data_c_n_rows_;
    for (unsigned int b = 0; b < solver->data_c_n_rows_; b += 8) {
      __builtin_prefetch(ptr + b, 0, 1);
    }
#endif
  }
};

}  // namespace fastcpd::families

#endif  // FASTCPD_FAMILIES_MGAUSSIAN_H_
