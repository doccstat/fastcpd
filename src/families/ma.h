#ifndef FASTCPD_FAMILIES_MA_H_
#define FASTCPD_FAMILIES_MA_H_

#include "fastcpd_family.h"

namespace fastcpd {
namespace families {

class MaFamily : public BaseFamily {
 public:
  static constexpr bool has_optimized_run = false;

  template <typename T>
  static void GetNllPeltValue(T* instance, unsigned int const segment_start,
                              unsigned int const segment_end, bool const cv,
                              std::optional<arma::colvec> const& start) {
    GetNllPelt(instance, segment_start, segment_end, cv, start);
  }

  template <typename T>
  static void GetNllPelt(T* instance, unsigned int const segment_start,
                         unsigned int const segment_end, bool const cv,
                         std::optional<arma::colvec> const& start) {
#ifndef NO_RCPP
    arma::mat const data_segment = instance->data_.rows(segment_start, segment_end);
    Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
    Rcpp::Function arima = stats["arima"];

    try {
      Rcpp::List out = arima(
          Rcpp::Named("x") = data_segment.col(0),
          Rcpp::Named("order") =
              Rcpp::NumericVector::create(0, 0, instance->order_(1)),
          Rcpp::Named("method") = "ML", Rcpp::Named("include.mean") = false);

      instance->result_coefficients_ = arma::zeros<arma::colvec>(instance->order_(1) + 1);
      instance->result_coefficients_.rows(0, instance->order_(1) - 1) =
          Rcpp::as<arma::colvec>(out["coef"]);
      instance->result_coefficients_(instance->order_(1)) = Rcpp::as<double>(out["sigma2"]);
      instance->result_residuals_ = arma::mat(Rcpp::as<arma::colvec>(out["residuals"]));
      instance->result_value_ = -Rcpp::as<double>(out["loglik"]);
    } catch (std::exception const& e) {
      Rcpp::warning("MA model fitting failed: %s", e.what());
      instance->result_coefficients_ = arma::zeros<arma::colvec>(instance->order_(1) + 1);
      instance->result_value_ = data_segment.n_rows * 10.0;
      instance->result_residuals_ = data_segment.col(0);
      instance->result_coefficients_(instance->order_(1)) = arma::var(data_segment.col(0));
    }
#endif
  }

  template <typename T>
  static double GetNllSen(T* instance, unsigned int const segment_start,
                          unsigned int const segment_end,
                          arma::colvec const& theta) {
    arma::mat data_segment = instance->data_.rows(segment_start, segment_end);
    unsigned int const q = instance->order_(1);
    arma::colvec reversed_theta = arma::reverse(theta);
    if (data_segment.n_rows < q + 1) {
      return 0;
    }
    arma::colvec variance_term = arma::zeros(data_segment.n_rows);
    for (unsigned int i = q; i < data_segment.n_rows; i++) {
      variance_term(i) =
          data_segment(i, 0) -
          arma::dot(reversed_theta.rows(1, q), variance_term.rows(i - q, i - 1));
    }
    return (std::log(2.0 * M_PI) + std::log(theta(q))) *
               (data_segment.n_rows - 2) / 2.0 +
           arma::dot(variance_term, variance_term) / 2.0 / theta(q);
  }

  template <typename T>
  static arma::colvec GetGradient(T* instance, unsigned int const segment_start,
                                  unsigned int const segment_end,
                                  arma::colvec const& theta) {
    arma::mat const data_segment = instance->data_.rows(segment_start, segment_end);
    unsigned int const segment_length = segment_end - segment_start + 1;
    unsigned int const q = instance->order_(1);
    arma::mat reversed_data = arma::reverse(data_segment, 0);
    arma::colvec reversed_theta = arma::reverse(theta);
    if (segment_length < q + 1) {
      return arma::ones(theta.n_elem);
    }
    arma::colvec variance_term = arma::zeros(segment_length);
    for (unsigned int i = q; i < segment_length; i++) {
      variance_term(i) =
          data_segment(i, 0) -
          arma::dot(reversed_theta.rows(1, q), variance_term.rows(i - q, i - 1));
    }
    arma::colvec reversed_variance_term = arma::reverse(variance_term);
    arma::mat psi_coefficient = arma::zeros(segment_length, q);
    for (unsigned int i = q; i < segment_length; i++) {
      psi_coefficient.row(i) =
          -reversed_variance_term
               .rows(segment_length - i, segment_length - i + q - 1)
               .t() -
          reversed_theta.rows(1, q).t() * psi_coefficient.rows(i - q, i - 1);
    }
    arma::colvec gradient = arma::zeros(q + 1);
    gradient.rows(0, q - 1) = psi_coefficient.row(segment_length - 1).t() *
                              variance_term(segment_length - 1) / theta(q);
    gradient(q) =
        1.0 / 2.0 / theta(q) -
        pow(variance_term(segment_length - 1), 2) / 2.0 / pow(theta(q), 2);
    return gradient;
  }

  template <typename T>
  static arma::mat GetHessian(T* instance, unsigned int const segment_start,
                              unsigned int const segment_end,
                              arma::colvec const& theta) {
    arma::mat const data_segment = instance->data_.rows(segment_start, segment_end);
    unsigned int const segment_length = segment_end - segment_start + 1;
    unsigned int const q = instance->order_(1);
    arma::mat reversed_data = arma::reverse(data_segment, 0);
    arma::colvec reversed_theta = arma::reverse(theta);
    if (segment_length < q + 1) {
      return arma::eye(theta.n_elem, theta.n_elem);
    }
    arma::colvec variance_term = arma::zeros(segment_length);
    for (unsigned int i = q; i < segment_length; i++) {
      variance_term(i) =
          data_segment(i, 0) -
          arma::dot(reversed_theta.rows(1, q), variance_term.rows(i - q, i - 1));
    }
    arma::colvec reversed_variance_term = arma::reverse(variance_term);
    arma::mat psi_coefficient = arma::zeros(segment_length, q);
    for (unsigned int i = q; i < segment_length; i++) {
      psi_coefficient.row(i) =
          -reversed_variance_term
               .rows(segment_length - i, segment_length - i + q - 1)
               .t() -
          reversed_theta.rows(1, q).t() * psi_coefficient.rows(i - q, i - 1);
    }
    arma::mat reversed_coef_psi = arma::reverse(psi_coefficient, 0);
    arma::cube psi_psi_coefficient = arma::zeros(q, q, segment_length);
    for (unsigned int i = q; i < segment_length; i++) {
      arma::mat psi_psi_coefficient_part = arma::zeros(q, q);
      for (unsigned int j = 1; j <= q; j++) {
        psi_psi_coefficient_part +=
            psi_psi_coefficient.slice(i - j) * theta(j - 1);
      }
      psi_psi_coefficient.slice(i) =
          -reversed_coef_psi.rows(segment_length - i,
                                  segment_length - i + q - 1) -
          reversed_coef_psi.rows(segment_length - i, segment_length - i + q - 1)
              .t() -
          psi_psi_coefficient_part;
    }
    arma::mat hessian = arma::zeros(q + 1, q + 1);
    hessian.submat(0, 0, q - 1, q - 1) =
        (psi_coefficient.row(segment_length - 1).t() *
             psi_coefficient.row(segment_length - 1) +
         psi_psi_coefficient.slice(segment_length - 1) *
             variance_term(segment_length - 1)) /
        theta(q);
    hessian.submat(0, q, q - 1, q) =
        -psi_coefficient.row(segment_length - 1).t() *
        variance_term(segment_length - 1) / theta(q) / theta(q);
    hessian.submat(q, 0, q, q - 1) = hessian.submat(0, q, q - 1, q).t();
    hessian(q, q) = pow(variance_term(segment_length - 1), 2) / pow(theta(q), 3) -
                    1.0 / 2.0 / pow(theta(q), 2);
    return hessian;
  }

  template <typename Solver>
  static void CreateSenParameters(Solver* solver) {
    solver->coefficients_.col(0) = solver->segment_coefficients_.row(0).t();
    solver->coefficients_sum_.col(0) = solver->segment_coefficients_.row(0).t();
    solver->hessian_.slice(0) =
        arma::zeros<arma::mat>(solver->parameters_count_, solver->parameters_count_);
  }

  template <typename Solver>
  static void UpdateSenParameters(Solver* solver) {
    int const segment_index =
        arma::index_max(arma::find(solver->segment_indices_ <= solver->t - 1));
    arma::colvec coef_add = solver->segment_coefficients_.row(segment_index).t();
    arma::mat hessian_new =
        MaFamily::GetHessian(solver, 0, solver->t - 1, coef_add);
    std::memcpy(solver->coefficients_.colptr(solver->pruned_set_size_ - 1),
                coef_add.memptr(),
                sizeof(double) * solver->parameters_count_);
    std::memcpy(solver->coefficients_sum_.colptr(solver->pruned_set_size_ - 1),
                coef_add.memptr(),
                sizeof(double) * solver->parameters_count_);
    std::memcpy(
        solver->hessian_.slice(solver->pruned_set_size_ - 1).memptr(),
        hessian_new.memptr(),
        sizeof(double) * solver->parameters_count_ * solver->parameters_count_);
  }
};

} // namespace families
} // namespace fastcpd

#endif
