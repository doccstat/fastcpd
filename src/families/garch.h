#ifndef FASTCPD_FAMILIES_GARCH_H_
#define FASTCPD_FAMILIES_GARCH_H_

#include "fastcpd_family.h"
#include "ref_tseries.h"

namespace fastcpd::families {

struct GarchFamily : BaseFamily {
  static constexpr const char* name = "garch";
  static constexpr bool has_optimized_run = false;
  static constexpr bool is_pelt_only = true;

  static unsigned int GetDataNDims(arma::mat const& data) {
    return data.n_cols;
  }

  template <typename Solver>
  static void GetNllPelt(Solver* solver, unsigned int const segment_start,
                         unsigned int const segment_end, bool const cv,
                         std::optional<arma::colvec> const& start) {
#ifndef NO_RCPP
    arma::colvec series = solver->data_.rows(segment_start, segment_end).col(0);
    Rcpp::List out = garch(series, solver->order_);
    solver->result_coefficients_ = Rcpp::as<arma::colvec>(out["coef"]);
    solver->result_residuals_ = arma::mat(Rcpp::as<arma::colvec>(out["residuals"]));
    solver->result_value_ = Rcpp::as<double>(out["n.likeli"]);
#endif
  }

  template <typename Solver>
  static void GetNllPeltValue(Solver* solver, unsigned int const segment_start,
                              unsigned int const segment_end, bool const cv,
                              std::optional<arma::colvec> const& start) {
    solver->result_value_ = GetNllPeltValueFast<-1>(solver, segment_start, segment_end);
  }

  template <int kNDims, typename Solver>
  static double GetNllPeltValueFast(Solver* solver, unsigned int const segment_start,
                                    unsigned int const segment_end) {
    static const std::optional<arma::colvec> empty_opt;
    GetNllPelt(solver, segment_start, segment_end, false, empty_opt);
    return solver->result_value_;
  }
};

}  // namespace fastcpd::families

#endif  // FASTCPD_FAMILIES_GARCH_H_
