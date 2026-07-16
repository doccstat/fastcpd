#ifndef FASTCPD_FAMILIES_ARIMA_H_
#define FASTCPD_FAMILIES_ARIMA_H_

#include <limits>
#include <utility>

#include "families/arma.h"
#include "fastcpd_family.h"

namespace fastcpd::families {

// ARIMA(p, d, q) change-point detection with segment-local differencing.
//
// Each candidate segment is differenced d times independently before the
// shared zero-mean ArmaFamily likelihood is evaluated. Differencing within
// the segment avoids introducing an artificial observation across a proposed
// change-point boundary and keeps returned change points in the coordinate
// system of the original, undifferenced series.
struct ArimaFamily : BaseFamily {
  static constexpr const char* name = "arima";
  static constexpr bool has_optimized_run = false;
  static constexpr bool is_pelt_only = true;

  template <typename Solver>
  static void GetNllPelt(Solver* solver, unsigned int const segment_start,
                         unsigned int const segment_end, bool const cv,
                         std::optional<arma::colvec> const& start) {
    FitSegment(solver, segment_start, segment_end, true);
  }

  template <typename Solver>
  static void GetNllPeltValue(
      Solver* solver, unsigned int const segment_start,
      unsigned int const segment_end, bool const cv,
      std::optional<arma::colvec> const& start) {
    solver->result_value_ =
        GetNllPeltValueFast<-1>(solver, segment_start, segment_end);
  }

  template <int kNDims, typename Solver>
  static double GetNllPeltValueFast(Solver* solver,
                                    unsigned int const segment_start,
                                    unsigned int const segment_end) {
    return FitSegment(solver, segment_start, segment_end, false);
  }

 private:
  struct ArmaFitAdapter {
    arma::mat data_;
    arma::colvec order_;
    arma::colvec result_coefficients_;
    arma::mat result_residuals_;
    double result_value_ = 0.0;

    ArmaFitAdapter(arma::colvec const& series, unsigned int const p,
                   unsigned int const q)
        : data_(series), order_(arma::colvec{static_cast<double>(p),
                                             static_cast<double>(q)}) {}
  };

  template <typename Solver>
  static double FitSegment(Solver* solver,
                           unsigned int const segment_start,
                           unsigned int const segment_end,
                           bool const keep_details) {
    unsigned int const p = static_cast<unsigned int>(solver->order_(0));
    unsigned int const d = static_cast<unsigned int>(solver->order_(1));
    unsigned int const q = static_cast<unsigned int>(solver->order_(2));
    unsigned int const segment_length = segment_end - segment_start + 1;

    if (segment_length <= d) {
      if (keep_details) {
        solver->result_coefficients_ =
            arma::zeros<arma::colvec>(solver->parameters_count_);
        solver->result_residuals_ = arma::mat(segment_length, 1);
        solver->result_residuals_.fill(
            std::numeric_limits<double>::quiet_NaN());
      }
      return solver->result_value_ = 0.0;
    }

    arma::colvec series =
        solver->data_.rows(segment_start, segment_end).col(0);
    arma::colvec const differenced = arma::diff(series, d);
    ArmaFitAdapter fit(differenced, p, q);
    static const std::optional<arma::colvec> empty_start;
    ArmaFamily::GetNllPelt(&fit, 0, differenced.n_elem - 1, false,
                           empty_start);

    if (keep_details) {
      solver->result_coefficients_ = std::move(fit.result_coefficients_);
      solver->result_residuals_ = arma::mat(segment_length, 1);
      solver->result_residuals_.fill(
          std::numeric_limits<double>::quiet_NaN());
      arma::uword const padding =
          segment_length - fit.result_residuals_.n_rows;
      if (fit.result_residuals_.n_rows > 0) {
        solver->result_residuals_.rows(padding, segment_length - 1) =
            fit.result_residuals_;
      }
    }
    return solver->result_value_ = fit.result_value_;
  }
};

}  // namespace fastcpd::families

#endif  // FASTCPD_FAMILIES_ARIMA_H_
