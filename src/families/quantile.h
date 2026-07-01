#ifndef FASTCPD_FAMILIES_QUANTILE_H_
#define FASTCPD_FAMILIES_QUANTILE_H_

#include <cmath>
#include "fastcpd_family.h"

namespace fastcpd::families {

struct QuantileFamily : BaseFamily {
  static constexpr const char* name = "quantile";
  static constexpr bool has_optimized_run = false;

  static unsigned int GetDataNDims(arma::mat const& data) {
    return data.n_cols - 1;
  }

  // IRLS quantile regression at quantile level tau.
  // Returns {coefficients, check-function cost = sum rho_tau(y - X*beta)}.
  static std::pair<arma::colvec, double> FitQuantileReg(
      arma::mat const& X, arma::colvec const& y, double const tau) {
    unsigned int const n = X.n_rows, p = X.n_cols;
    double const kEps = 1e-6;
    arma::colvec beta =
        arma::solve(X.t() * X + kEps * arma::eye(p, p), X.t() * y);
    for (unsigned int iter = 0; iter < 100; ++iter) {
      arma::colvec const r = y - X * beta;
      arma::colvec w(n);
      for (unsigned int i = 0; i < n; ++i) {
        w(i) = (r(i) > 0.0 ? tau : (1.0 - tau)) /
               std::max(std::abs(r(i)), kEps);
      }
      arma::mat WX = X;
      WX.each_col() %= w;
      arma::colvec beta_new =
          arma::solve(X.t() * WX + kEps * arma::eye(p, p), WX.t() * y);
      double const delta = arma::norm(beta_new - beta);
      beta = beta_new;
      if (delta < 1e-8 * (1.0 + arma::norm(beta))) break;
    }
    arma::colvec const r = y - X * beta;
    double cost = 0.0;
    for (unsigned int i = 0; i < n; ++i) {
      cost += r(i) >= 0.0 ? tau * r(i) : (tau - 1.0) * r(i);
    }
    return {beta, cost};
  }

  template <typename Solver>
  static void GetNllPelt(Solver* solver, unsigned int const segment_start,
                         unsigned int const segment_end, bool const cv,
                         std::optional<arma::colvec> const& start) {
    arma::mat const data_segment =
        solver->data_.rows(segment_start, segment_end);
    arma::colvec const y = data_segment.col(0);
    arma::mat const X = data_segment.cols(1, data_segment.n_cols - 1);
    double const tau = solver->order_(0);
    auto const [beta, cost] = FitQuantileReg(X, y, tau);
    solver->result_coefficients_ = beta;
    solver->result_residuals_ = arma::mat(y - X * beta);
    solver->result_value_ = cost;
  }

  template <typename Solver>
  static void GetNllPeltValue(Solver* solver, unsigned int const segment_start,
                              unsigned int const segment_end, bool const cv,
                              std::optional<arma::colvec> const& start) {
    GetNllPelt(solver, segment_start, segment_end, cv, start);
  }

  template <typename Solver>
  static double GetNllSen(Solver* solver, unsigned int const segment_start,
                          unsigned int const segment_end,
                          arma::colvec const& theta) {
    double const tau = solver->order_(0);
    double cost = 0.0;
    for (unsigned int row = segment_start; row <= segment_end; row++) {
      double const residual =
          solver->data_(row, 0) - DotTail(solver, row, theta);
      cost += residual >= 0.0 ? tau * residual : (tau - 1.0) * residual;
    }
    return cost;
  }

  template <typename Solver>
  static arma::colvec GetGradient(Solver* solver,
                                  unsigned int const segment_start,
                                  unsigned int const segment_end,
                                  arma::colvec const& theta) {
    double const y = solver->data_(segment_end, 0);
    double const tau = solver->order_(0);
    double const residual = y - DotTail(solver, segment_end, theta);
    double const indicator = (residual < 0.0) ? 1.0 : 0.0;
    return ScaledTail(solver, segment_end, indicator - tau);
  }

  template <typename Solver>
  static arma::mat GetHessian(Solver* solver, unsigned int const segment_start,
                              unsigned int const segment_end,
                              arma::colvec const& theta) {
    return TailOuter(solver, segment_end);
  }

  template <typename Solver>
  static void CreateSenParameters(Solver* solver) {
    solver->coefficients_.col(0) = solver->segment_coefficients_.row(0).t();
    solver->coefficients_sum_.col(0) = solver->segment_coefficients_.row(0).t();
    solver->hessian_.slice(0) =
        TailOuter(solver, 0) +
        solver->epsilon_in_hessian_ *
            arma::eye<arma::mat>(solver->parameters_count_,
                                 solver->parameters_count_);
  }

  template <typename Solver>
  static void UpdateSenParameters(Solver* solver) {
    unsigned int const segment_index = solver->segment_index_;
    arma::colvec const coef_add =
        solver->segment_coefficients_.row(segment_index).t();
    arma::mat const hessian_new =
        TailOuter(solver, solver->t - 1) +
        solver->epsilon_in_hessian_ *
            arma::eye<arma::mat>(solver->parameters_count_,
                                 solver->parameters_count_);
    std::memcpy(solver->coefficients_.colptr(solver->pruned_set_size_ - 1),
                coef_add.memptr(),
                sizeof(double) * solver->parameters_count_);
    std::memcpy(solver->coefficients_sum_.colptr(solver->pruned_set_size_ - 1),
                coef_add.memptr(),
                sizeof(double) * solver->parameters_count_);
    std::memcpy(solver->hessian_.slice(solver->pruned_set_size_ - 1).memptr(),
                hessian_new.memptr(),
                sizeof(double) * solver->parameters_count_ *
                    solver->parameters_count_);
  }
};

}  // namespace fastcpd::families

#endif  // FASTCPD_FAMILIES_QUANTILE_H_
