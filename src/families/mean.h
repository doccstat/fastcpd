#ifndef FASTCPD_FAMILIES_MEAN_H_
#define FASTCPD_FAMILIES_MEAN_H_

#include "fastcpd_family.h"

namespace fastcpd::families {

struct MeanFamily : BaseFamily {
  static constexpr const char* name = "mean";
  // Optimisation lives in GetNllPeltValue (O(1) cumulative-sum lookup).
  // The standard PELT loop is used; no separate Run() override is needed.
  static constexpr bool has_optimized_run = false;
  static constexpr bool is_pelt_only = true;

  // Data layout: TRANSPOSED relative to other families.
  //
  // Other families store data_c_ as (n+1) × (p+1) column-major, so element
  // [t, d] is at offset  t + d * (n+1).  Accessing all p+1 dimensions for a
  // given candidate s therefore jumps by n+1 doubles (≈ 8 GB for n = 10^9),
  // causing one TLB miss + DRAM miss per dimension per candidate.
  //
  // We instead return data_c_.t() — shape (p+1) × (n+1) — so element [d, t]
  // is at offset  d + t * (p+1).  All p+1 values for candidate s now live in
  // (p+1)*8 consecutive bytes: for p=1 that is one 16-byte SSE slot; for p=7
  // it fits in a single 64-byte cache line.  Hardware prefetch automatically
  // brings in the fixed endpoint (t moves forward by p+1 doubles per step),
  // and PrefetchCandidate hides latency for random-s reads.
  //
  // Consequence: data_c_n_rows_ == p+1 after the transpose.  The stride
  // between consecutive time indices is p+1 (tiny), not n+1 (huge).
  static arma::mat CreateDataC(arma::mat const& data,
                               arma::mat const& variance_estimate,
                               unsigned int const p_response = 0) {
    arma::mat data_c = data * arma::chol(arma::inv(variance_estimate)).t();
    data_c =
        arma::cumsum(arma::join_rows(data_c, arma::sum(arma::square(data_c), 1)));
    data_c = arma::join_cols(arma::zeros<arma::rowvec>(data_c.n_cols), data_c);
    return data_c.t();  // (p+1) × (n+1): time is the column index
  }

  static unsigned int GetDataNDims(arma::mat const& data) {
    return data.n_cols;
  }

  template <typename Solver>
  static void GetNllPelt(Solver* solver, unsigned int const segment_start,
                         unsigned int const segment_end, bool const cv,
                         std::optional<arma::colvec> const& start) {
    arma::mat data_segment = solver->data_.rows(segment_start, segment_end);
    solver->result_coefficients_ = arma::mean(data_segment, 0).t();
    solver->result_residuals_ =
        data_segment.each_row() - solver->result_coefficients_.t();
    GetNllPeltValue(solver, segment_start, segment_end, cv, start);
  }

  template <typename Solver>
  static void GetNllPeltValue(Solver* solver, unsigned int const segment_start,
                              unsigned int const segment_end, bool const cv,
                              std::optional<arma::colvec> const& start) {
    solver->result_value_ = GetNllPeltValueFast<-1>(solver, segment_start, segment_end);
  }

  // Core NLL computation.
  //
  // With the transposed layout, stride == data_c_n_rows_ == p+1.  For a given
  // segment [start, end]:
  //   end_ptr   = data_c + (end+1) * stride  — fixed for all candidates at t
  //   start_ptr = data_c + start   * stride  — varies per candidate
  // The inner loop end_ptr[d] - start_ptr[d] is stride-1 sequential access,
  // auto-vectorized to SSE2/AVX2/NEON by the compiler.
  //
  // When kNDims == 1 (p == 1), the loop degenerates to a single scalar
  // computation — emitted as two loads, one sub, one mul, no loop overhead.
  //
  // Alignment: Armadillo guarantees >= 16-byte alignment for all heap
  // allocations via posix_memalign / _aligned_malloc (matching its own
  // memory::mark_as_aligned which always hints 16).  We use 16 here for the
  // same reason: the hint is on the BASE pointer; end_ptr and start_ptr are
  // derived with a runtime offset so the compiler cannot propagate a stronger
  // alignment to them regardless.  No runtime branch — always 16.
  template <int kNDims, typename Solver>
  static double GetNllPeltValueFast(Solver* solver, unsigned int const segment_start,
                                    unsigned int const segment_end) {
#if defined(__GNUC__) || defined(__clang__)
    double const* __restrict__ data_c = static_cast<double const*>(
        __builtin_assume_aligned(solver->data_c_ptr_, 16));
#else
    double const* __restrict__ data_c = solver->data_c_ptr_;
#endif
    unsigned int const stride = solver->data_c_n_rows_;  // == p+1
    double const* const end_ptr =
        data_c + static_cast<std::size_t>(segment_end + 1) * stride;
    double const* const start_ptr =
        data_c + static_cast<std::size_t>(segment_start) * stride;
    unsigned int const segment_length = segment_end - segment_start + 1;
    if constexpr (kNDims == 1) {
      // p == 1: single scalar computation, zero loop overhead.
      // stride == p+1 == 2: end_ptr[0] = cumsum_x, end_ptr[1] = cumsum_x².
      double const diff0 = end_ptr[0] - start_ptr[0];
      return (end_ptr[1] - start_ptr[1] - diff0 * diff0 / segment_length) / 2.0;
    }
    double two_norm = 0;
    // No loop-carried memory dependency — end_ptr and start_ptr are distinct
    // non-aliased regions; ivdep lets the compiler emit fully unrolled SIMD.
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
    for (unsigned int d = 0; d < solver->parameters_count_; d++) {
      double const diff = end_ptr[d] - start_ptr[d];
      two_norm += diff * diff;
    }
    return (end_ptr[solver->parameters_count_] -
            start_ptr[solver->parameters_count_] -
            two_norm / segment_length) /
           2.0;
  }

  // Prefetch all p+1 values for candidate s before they are needed.
  // With the transposed layout they begin at data_c_ptr_ + s*(p+1).
  // Two unconditional prefetches — no runtime branch:
  //   line 0 covers dimensions 0..7   (64 bytes, p <= 7 fits entirely)
  //   line 1 covers dimensions 8..15  (second 64-byte cache line, p 8..15)
  // For p <= 7 the second prefetch reaches into the next candidate's first
  // line, warming it up slightly early — harmless and often beneficial.
  // For p > 15 add further lines here as needed.
  template <typename Solver>
  static void PrefetchCandidate(Solver* solver, unsigned int const s) {
#if defined(__GNUC__) || defined(__clang__)
    double const* const ptr =
        solver->data_c_ptr_ + static_cast<std::size_t>(s) * solver->data_c_n_rows_;
    __builtin_prefetch(ptr,     0, 1);
    __builtin_prefetch(ptr + 8, 0, 1);
#endif
  }
};

}  // namespace fastcpd::families

#endif  // FASTCPD_FAMILIES_MEAN_H_
