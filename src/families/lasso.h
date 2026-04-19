#ifndef FASTCPD_FAMILIES_LASSO_H_
#define FASTCPD_FAMILIES_LASSO_H_

#include "fastcpd_family.h"

namespace fastcpd::families {

struct LassoFamily : BaseFamily {
  static constexpr const char* name = "lasso";
  static constexpr bool has_optimized_run = false;
  static constexpr bool has_error_std_dev = true;
  static constexpr bool is_lasso = true;

  static unsigned int GetDataNDims(arma::mat const& data) {
    return data.n_cols - 1;
  }

  template <typename Solver>
  static void CreateSenParameters(Solver* solver) {
    solver->coefficients_.col(0) = solver->segment_coefficients_.row(0).t();
    solver->coefficients_sum_.col(0) = solver->segment_coefficients_.row(0).t();
    solver->hessian_.slice(0) =
        solver->data_.row(0).tail(solver->parameters_count_).t() *
            solver->data_.row(0).tail(solver->parameters_count_) +
        solver->epsilon_in_hessian_ *
            arma::eye<arma::mat>(solver->parameters_count_,
                                 solver->parameters_count_);
  }

  template <typename Solver>
  static void UpdateSenParameters(Solver* solver) {
    int const segment_index =
        arma::index_max(arma::find(solver->segment_indices_ <= solver->t - 1));
    arma::colvec coef_add = solver->segment_coefficients_.row(segment_index).t();
    arma::rowvec const x =
        solver->data_.row(solver->t - 1).tail(solver->parameters_count_);
    arma::mat hessian_new =
        x.t() * x + solver->epsilon_in_hessian_ *
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
#ifdef NO_RCPP
    // Pure-C++ fallback (Python binding): use OLS instead of glmnet.
    // This gives correct initialization for the SEN gradient-descent path
    // even though the regularisation is not applied here.
    arma::mat const data_segment =
        solver->data_.rows(segment_start, segment_end);
    arma::mat const x = data_segment.cols(1, data_segment.n_cols - 1);
    arma::colvec const y = data_segment.col(0);
    unsigned int const p = x.n_cols;
    if (segment_start == segment_end || x.n_rows <= p) {
      solver->result_coefficients_ = arma::zeros<arma::colvec>(p);
      solver->result_residuals_ = y;
      solver->result_value_ = arma::dot(y, y) / 2.0;
    } else {
      // OLS: β = (X'X + ε·I)⁻¹ X'y  (ridge for numerical stability)
      arma::mat const XtX =
          x.t() * x + solver->epsilon_in_hessian_ * arma::eye(p, p);
      solver->result_coefficients_ = arma::solve(XtX, x.t() * y);
      solver->result_residuals_ = y - x * solver->result_coefficients_;
      solver->result_value_ =
          arma::dot(solver->result_residuals_, solver->result_residuals_) / 2.0;
    }
#else
    if (segment_start == segment_end) {
      solver->result_coefficients_ =
          arma::zeros<arma::colvec>(solver->data_.n_cols - 1);
      solver->result_residuals_ = arma::zeros<arma::mat>(1, 1);
      solver->result_value_ = 0;
    } else if (cv) {
      arma::mat const data_segment =
          solver->data_.rows(segment_start, segment_end);
      Rcpp::Environment glmnet = Rcpp::Environment::namespace_env("glmnet"),
                        stats = Rcpp::Environment::namespace_env("stats");
      Rcpp::Function cv_glmnet = glmnet["cv.glmnet"],
                     predict_glmnet = glmnet["predict.glmnet"],
                     deviance = stats["deviance"];
      Rcpp::List out =
          cv_glmnet(data_segment.cols(1, data_segment.n_cols - 1),
                    data_segment.col(0), Rcpp::Named("family") = "gaussian");
      arma::colvec index_vec = Rcpp::as<arma::colvec>(out["index"]),
                   values = Rcpp::as<arma::colvec>(deviance(out["glmnet.fit"]));
      Rcpp::S4 out_coef = predict_glmnet(
          out["glmnet.fit"], Rcpp::Named("s") = out["lambda.1se"],
          Rcpp::Named("type") = "coefficients", Rcpp::Named("exact") = false);
      arma::colvec glmnet_i = Rcpp::as<arma::colvec>(out_coef.slot("i"));
      arma::colvec glmnet_x = Rcpp::as<arma::colvec>(out_coef.slot("x"));
      solver->result_coefficients_ =
          arma::zeros<arma::colvec>(data_segment.n_cols - 1);
      for (unsigned int i = 1; i < glmnet_i.n_elem; i++) {
        solver->result_coefficients_(glmnet_i(i) - 1) = glmnet_x(i);
      }
      solver->result_residuals_ = arma::mat();
      solver->result_value_ = values(index_vec(1) - 1);
    } else {
      arma::mat const data_segment =
          solver->data_.rows(segment_start, segment_end);
      Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats"),
                        glmnet = Rcpp::Environment::namespace_env("glmnet");
      Rcpp::Function deviance = stats["deviance"], glmnet_ = glmnet["glmnet"],
                     predict_glmnet = glmnet["predict.glmnet"];
      Rcpp::List out = glmnet_(
          data_segment.cols(1, data_segment.n_cols - 1), data_segment.col(0),
          Rcpp::Named("family") = "gaussian",
          Rcpp::Named("lambda") = solver->lasso_penalty_base_ /
                                  std::sqrt(segment_end - segment_start + 1));
      Rcpp::S4 out_par = out["beta"];
      arma::colvec par_i = Rcpp::as<arma::colvec>(out_par.slot("i"));
      arma::colvec par_x = Rcpp::as<arma::colvec>(out_par.slot("x"));
      solver->result_coefficients_ =
          arma::zeros<arma::colvec>(data_segment.n_cols - 1);
      for (unsigned int i = 0; i < par_i.n_elem; i++) {
        solver->result_coefficients_(par_i(i)) = par_x(i);
      }
      double value = Rcpp::as<double>(deviance(out));
      arma::colvec fitted_values = Rcpp::as<arma::colvec>(predict_glmnet(
          out, data_segment.cols(1, data_segment.n_cols - 1),
          Rcpp::Named("s") = solver->lasso_penalty_base_ /
                             std::sqrt(segment_end - segment_start + 1)));
      solver->result_residuals_ = data_segment.col(0) - fitted_values;
      solver->result_value_ = value / 2;
    }
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
    return arma::accu(arma::square(y - x * theta)) / 2 +
           solver->lasso_penalty_base_ /
               std::sqrt(segment_end - segment_start + 1) *
               arma::accu(arma::abs(theta));
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
    return -(y - arma::as_scalar(x * theta)) * x.t();
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
    return x.t() * x;
  }
};

}  // namespace fastcpd::families

#endif  // FASTCPD_FAMILIES_LASSO_H_
