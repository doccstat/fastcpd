#ifndef FASTCPD_FAMILIES_GAUSSIAN_H_
#define FASTCPD_FAMILIES_GAUSSIAN_H_

#include "fastcpd_family.h"
#include "fastcpd_glm.h"

namespace fastcpd::families {

struct GaussianFamily : BaseFamily {
  static constexpr const char* name = "gaussian";
  static constexpr bool has_optimized_run = false;
  static constexpr bool has_error_std_dev = true;

  static unsigned int GetDataNDims(arma::mat const& data) {
    return data.n_cols - 1;
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
    arma::colvec coef_add = solver->segment_coefficients_.row(segment_index).t();
    arma::mat hessian_new =
        TailOuter(solver, solver->t - 1) +
        solver->epsilon_in_hessian_ *
            arma::eye<arma::mat>(solver->parameters_count_,
                                 solver->parameters_count_);

    std::memcpy(solver->coefficients_.colptr(solver->pruned_set_size_ - 1),
                coef_add.memptr(), sizeof(double) * solver->parameters_count_);
    std::memcpy(solver->coefficients_sum_.colptr(solver->pruned_set_size_ - 1),
                coef_add.memptr(), sizeof(double) * solver->parameters_count_);
    std::memcpy(solver->hessian_.slice(solver->pruned_set_size_ - 1).memptr(),
                hessian_new.memptr(),
                sizeof(double) * solver->parameters_count_ *
                    solver->parameters_count_);
  }

  template <typename Solver>
  static void GetNllPelt(Solver* solver, unsigned int const segment_start,
                         unsigned int const segment_end, bool const cv,
                         std::optional<arma::colvec> const& start) {
    arma::mat const data_segment = solver->data_.rows(segment_start, segment_end);
    arma::colvec y = data_segment.col(0);
    arma::mat x = data_segment.cols(1, data_segment.n_cols - 1);
    fastcpd::glm::FitGlmResult out;
    if (!start.has_value()) {
      out = fastcpd::glm::FitGlm(x, y, name);
    } else {
      out = fastcpd::glm::FitGlm(x, y, name, start);
    }
    solver->result_coefficients_ = out.coefficients;
    solver->result_residuals_ = arma::mat(out.residuals);
    solver->result_value_ = out.deviance / 2;
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
    double cost = 0.0;
    for (unsigned int row = segment_start; row <= segment_end; row++) {
      double const residual = solver->data_(row, 0) - DotTail(solver, row, theta);
      cost += residual * residual;
    }
    return cost / 2.0;
  }

  template <typename Solver>
  static arma::colvec GetGradient(Solver* solver,
                                  unsigned int const segment_start,
                                  unsigned int const segment_end,
                                  arma::colvec const& theta) {
    double const y = solver->data_(segment_end, 0);
    double const eta = DotTail(solver, segment_end, theta);
    return ScaledTail(solver, segment_end, eta - y);
  }

  template <typename Solver>
  static arma::mat GetHessian(Solver* solver, unsigned int const segment_start,
                              unsigned int const segment_end,
                              arma::colvec const& theta) {
    return TailOuter(solver, segment_end);
  }
};

}  // namespace fastcpd::families

#endif  // FASTCPD_FAMILIES_GAUSSIAN_H_
