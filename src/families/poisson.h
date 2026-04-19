#ifndef FASTCPD_FAMILIES_POISSON_H_
#define FASTCPD_FAMILIES_POISSON_H_

#include "fastcpd_family.h"
#ifndef NO_RCPP
#include "ref_fastglm_fit_glm.h"
#endif

namespace fastcpd::families {

struct PoissonFamily : BaseFamily {
  static constexpr const char* name = "poisson";
  static constexpr bool has_optimized_run = false;
  static constexpr bool supports_warm_start = true;

  static unsigned int GetDataNDims(arma::mat const& data) {
    return data.n_cols - 1;
  }

  template <typename Solver>
  static void CreateSenParameters(Solver* solver) {
    solver->coefficients_.col(0) = solver->segment_coefficients_.row(0).t();
    solver->coefficients_sum_.col(0) = solver->segment_coefficients_.row(0).t();
    solver->hessian_.slice(0) =
        (solver->data_.row(0).tail(solver->parameters_count_).t() *
         solver->data_.row(0).tail(solver->parameters_count_)) *
        std::exp(arma::dot(
            solver->coefficients_.col(0),
            solver->data_.row(0).tail(solver->parameters_count_)));
  }

  template <typename Solver>
  static void UpdateSenParameters(Solver* solver) {
    int const segment_index =
        arma::index_max(arma::find(solver->segment_indices_ <= solver->t - 1));
    arma::colvec coef_add = solver->segment_coefficients_.row(segment_index).t();
    arma::rowvec const x =
        solver->data_.row(solver->t - 1).tail(solver->parameters_count_);
    arma::mat hessian_new = (x.t() * x) * std::exp(arma::dot(coef_add, x));

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
#ifndef NO_RCPP
    arma::mat const data_segment = solver->data_.rows(segment_start, segment_end);
    arma::colvec y = data_segment.col(0);
    arma::mat x = data_segment.cols(1, data_segment.n_cols - 1);
    Rcpp::List out;
    if (!start.has_value()) {
      out = fastglm(x, y, name);
    } else {
      out = fastglm(x, y, name, start);
    }
    solver->result_coefficients_ = Rcpp::as<arma::colvec>(out["coefficients"]);
    solver->result_residuals_ =
        arma::mat(Rcpp::as<arma::colvec>(out["residuals"]));
    solver->result_value_ = Rcpp::as<double>(out["deviance"]) / 2;
#endif
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
    arma::mat data_segment = solver->data_.rows(segment_start, segment_end);
    arma::colvec y = data_segment.col(0);
    arma::mat x = data_segment.cols(1, data_segment.n_cols - 1);
    arma::colvec u = x * theta;
    arma::colvec y_factorial(y.n_elem);
    for (unsigned int i = 0; i < y.n_elem; i++) {
      double log_factorial = 0;
      for (int j = 1; j <= y(i); ++j) {
        log_factorial += std::log(j);
      }
      y_factorial(i) = log_factorial;
    }
    return arma::accu(-y % u + arma::exp(u) + y_factorial);
  }

  template <typename Solver>
  static arma::colvec GetGradient(Solver* solver,
                                  unsigned int const segment_start,
                                  unsigned int const segment_end,
                                  arma::colvec const& theta) {
    arma::mat const data_segment =
        solver->data_.rows(segment_start, segment_end);
    unsigned int const segment_length = segment_end - segment_start + 1;
    arma::rowvec new_data = data_segment.row(segment_length - 1);
    arma::rowvec x = new_data.tail(new_data.n_elem - 1);
    double y = new_data(0);
    return -(y - std::exp(arma::as_scalar(x * theta))) * x.t();
  }

  template <typename Solver>
  static arma::mat GetHessian(Solver* solver, unsigned int const segment_start,
                              unsigned int const segment_end,
                              arma::colvec const& theta) {
    arma::mat const data_segment =
        solver->data_.rows(segment_start, segment_end);
    unsigned int const segment_length = segment_end - segment_start + 1;
    arma::rowvec new_data = data_segment.row(segment_length - 1);
    arma::rowvec x = new_data.tail(new_data.n_elem - 1);
    double prob = std::exp(arma::as_scalar(x * theta));
    return (x.t() * x) * std::min(arma::as_scalar(prob), 1e10);
  }
};

}  // namespace fastcpd::families

#endif  // FASTCPD_FAMILIES_POISSON_H_
