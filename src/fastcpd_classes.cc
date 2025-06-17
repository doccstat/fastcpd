#include "fastcpd.h"
#ifndef NO_RCPP
#include "ref_fastglm_fit_glm.h"
#include "ref_tseries.h"
#endif
#include <cmath>
#include <optional>

using ::arma::abs;
using ::arma::accu;
using ::arma::as_scalar;
using ::arma::colvec;
using ::arma::cube;
using ::arma::diff;
using ::arma::dot;
using ::arma::eye;
using ::arma::find;
using ::arma::floor;
using ::arma::index_max;
using ::arma::index_min;
using ::arma::join_cols;
using ::arma::join_rows;
using ::arma::join_slices;
using ::arma::linspace;
using ::arma::log_det_sympd;
using ::arma::mat;
using ::arma::max;
using ::arma::mean;
using ::arma::min;
using ::arma::norm;
using ::arma::ones;
using ::arma::reverse;
using ::arma::rowvec;
using ::arma::sign;
using ::arma::solve;
using ::arma::sort;
using ::arma::square;
using ::arma::trace;
using ::arma::ucolvec;
using ::arma::unique;
using ::arma::vec;
using ::arma::zeros;
using ::std::make_unique;
using ::std::memcpy;
using ::std::pow;
using ::std::string;
using ::std::string_view;
using ::std::unique_ptr;
using ::std::unordered_map;
using ::std::vector;

namespace fastcpd::classes {

Fastcpd::Fastcpd(
    double const beta,
#ifndef NO_RCPP
    std::optional<Rcpp::Function> const& cost,
#endif
    std::function<double(arma::mat)> const& cost_pelt,
    std::function<double(arma::mat, arma::colvec)> const& cost_sen,
    std::string const& cost_adjustment,
    std::function<arma::colvec(arma::mat, arma::colvec)> const& cost_gradient,
    std::function<arma::mat(arma::mat, arma::colvec)> const& cost_hessian,
    bool const cp_only, arma::mat const& data, double const epsilon,
    std::string const& family,
    std::function<unsigned int(unsigned int)> const& multiple_epochs_function,
    arma::colvec const& line_search, arma::colvec const& lower,
    double const momentum_coef, arma::colvec const& order, int const p,
    unsigned int const p_response, double const pruning_coef,
    bool const r_progress, int const segment_count, double const trim,
    arma::colvec const& upper, double const vanilla_percentage,
    arma::mat const& variance_estimate, bool const warm_start)
    : active_coefficients_count_(colvec(segment_count)),
      beta_(beta),
      change_points_(zeros<colvec>(data.n_rows + 1)),
      coefficients_(mat(p, data.n_rows + 1)),
      coefficients_sum_(mat(p, data.n_rows + 1)),
      cost_adjustment_(cost_adjustment),
#ifndef NO_RCPP
      cost_function_([&]() -> std::optional<Rcpp::Function> {
        if (family == "custom") {
          return cost;
        }
        return std::nullopt;
      }()),
#endif
      cost_function_pelt_(cost_pelt),
      cost_function_sen_(cost_sen),
      cost_gradient_(cost_gradient),
      cost_hessian_(cost_hessian),
      cp_only_(cp_only),
      data_(data),
      data_c_([&]() -> mat {
        if (family == "mean") {
          mat data_c = data * chol(inv(variance_estimate)).t();
          data_c = cumsum(join_rows(data_c, sum(square(data_c), 1)));
          data_c = join_cols(zeros<rowvec>(data_c.n_cols), data_c);
          return data_c;
        } else if (family == "variance") {
          mat data_c = data;
          data_c.each_row() -= mean(data_c, 0);
          mat data_crossprod(data.n_cols * data.n_cols, data.n_rows);
          for (unsigned int i = 0; i < data.n_rows; i++) {
            data_crossprod.col(i) =
                vectorise(data_c.row(i).t() * data_c.row(i));
          }
          data_c = cumsum(data_crossprod.t());
          data_c = join_cols(zeros<rowvec>(data_c.n_cols), data_c);
          return data_c;
        } else if (family == "meanvariance") {
          mat data_crossprod(data.n_cols * data.n_cols, data.n_rows);
          for (unsigned int i = 0; i < data.n_rows; i++) {
            data_crossprod.col(i) = vectorise(data.row(i).t() * data.row(i));
          }
          mat data_c = cumsum(join_rows(data, data_crossprod.t()));
          data_c = join_cols(zeros<rowvec>(data_c.n_cols), data_c);
          return data_c;
        }
        return data;
      }()),
      data_c_n_cols_(data_c_.n_cols),
      data_c_n_rows_(data_c_.n_rows),
      data_c_ptr_(data_c_.memptr()),
      data_n_dims_([&]() -> unsigned int {
        if (family == "mean" || family == "variance" ||
            family == "meanvariance" || family == "ar" || family == "arma" ||
            family == "arima" || family == "garch" || family == "var" ||
            family == "ma") {
          return data.n_cols;
        }
        return data.n_cols - 1;
      }()),
      data_n_cols_(data_.n_cols),
      data_n_rows_(data_.n_rows),
      epsilon_in_hessian_(epsilon),
      error_standard_deviation_(colvec(segment_count)),
      family_(family),
      get_gradient_([&]() -> arma::colvec (Fastcpd::*)(unsigned int const,
                                                       unsigned int const,
                                                       arma::colvec const&) {
        auto it = family_function_map_.find(family);
        if (it != family_function_map_.end()) {
          FunctionSet const& func_set = it->second;
          return func_set.gradient;
        }
        return &Fastcpd::GetGradientCustom;
      }()),
      get_hessian_([&]() -> arma::mat (Fastcpd::*)(unsigned int const,
                                                   unsigned int const,
                                                   arma::colvec const&) {
        auto it = family_function_map_.find(family);
        if (it != family_function_map_.end()) {
          FunctionSet const& func_set = it->second;
          return func_set.hessian;
        }
        return &Fastcpd::GetHessianCustom;
      }()),
      get_nll_pelt_([&]() -> void (Fastcpd::*)(
                              unsigned int const, unsigned int const,
                              bool const, std::optional<arma::colvec> const&) {
        auto it = family_function_map_.find(family);
        if (it != family_function_map_.end()) {
          FunctionSet const& func_set = it->second;
          return func_set.nll_pelt;
        }
        return &Fastcpd::GetNllPeltCustom;
      }()),
      get_nll_pelt_value_(
          [&]() -> void (Fastcpd::*)(unsigned int const, unsigned int const,
                                     bool const,
                                     std::optional<arma::colvec> const&) {
            auto it = family_function_map_.find(family);
            if (it != family_function_map_.end()) {
              FunctionSet const& func_set = it->second;
              return func_set.nll_pelt_value;
            }
            return &Fastcpd::GetNllPeltCustom;
          }()),
      get_nll_sen_([&]() -> double (Fastcpd::*)(unsigned int const,
                                                unsigned int const,
                                                arma::colvec const&) {
        auto it = family_function_map_.find(family);
        if (it != family_function_map_.end()) {
          FunctionSet const& func_set = it->second;
          return func_set.nll_sen;
        }
        return &Fastcpd::GetNllSenCustom;
      }()),
      hessian_(cube(p, p, data.n_rows + 1)),
      line_search_(line_search),
      momentum_(colvec(p)),
      momentum_coef_(momentum_coef),
      multiple_epochs_function_(multiple_epochs_function),
      objective_function_values_(colvec(data_n_rows_ + 1)),
      objective_function_values_candidates_(colvec(data_n_rows_ + 1)),
      objective_function_values_candidates_ptr_(
          objective_function_values_candidates_.memptr()),
      order_(order),
      parameters_count_(p),
      parameters_lower_bound_(lower),
      parameters_upper_bound_(upper),
      pruned_left_(ucolvec(data_n_rows_ + 1)),
      pruned_set_(zeros<ucolvec>(data_n_rows_ + 1)),
      pruning_coefficient_(pruning_coef),
      r_progress_(r_progress),
      regression_response_count_(p_response),
#ifndef NO_RCPP
      rProgress_(make_unique<RProgress::RProgress>(kRProgress, data_n_rows_)),
#endif
      segment_coefficients_(mat(segment_count, p)),
      segment_count_(segment_count),
      segment_indices_(round(linspace(0, data_n_rows_, segment_count + 1))),
      trim_(trim),
      use_warm_start_(warm_start),
      vanilla_percentage_(vanilla_percentage),
      variance_estimate_(variance_estimate),
      warm_start_(zeros<mat>(p, data_n_rows_)) {
  // TODO(doccstat): Store environment functions from R.
}

std::tuple<arma::colvec, arma::colvec, arma::colvec, arma::mat, arma::mat>
Fastcpd::Run() {
  CreateSegmentStatistics();
  CreateSenParameters();
  pruned_set_(1) = 1;
  objective_function_values_.fill(arma::datum::inf);
  objective_function_values_(0) = -beta_;
  objective_function_values_(1) = GetCostValue(0, 0);

  // TODO(doccstat): Investigate if the following branches can be merged into
  // `fastcpd_class_nll.cc`.
  if (family_ == "mean" && cost_adjustment_ == "MBIC") {
    double two_norm;
    unsigned int i, pi;

    for (t = 2; t <= data_n_rows_; t++) {
      for (i = 0; i < pruned_set_size_; i++) {
        two_norm = (data_c_ptr_[t] - data_c_ptr_[pruned_set_[i]]) *
                   (data_c_ptr_[t] - data_c_ptr_[pruned_set_[i]]);
        for (pi = 1; pi < parameters_count_; pi++) {
          two_norm += (data_c_ptr_[t + data_c_n_rows_ * pi] -
                       data_c_ptr_[pruned_set_[i] + data_c_n_rows_ * pi]) *
                      (data_c_ptr_[t + data_c_n_rows_ * pi] -
                       data_c_ptr_[pruned_set_[i] + data_c_n_rows_ * pi]);
        }
        objective_function_values_candidates_ptr_[i] =
            objective_function_values_[pruned_set_[i]] +
            ((data_c_ptr_[t + data_c_n_rows_ * parameters_count_] -
              data_c_ptr_[pruned_set_[i] +
                          data_c_n_rows_ * parameters_count_]) -
             two_norm / (t - pruned_set_[i])) /
                2.0 +
            std::log(t - pruned_set_[i]) / 2.0 + beta_;
      }

      objective_function_values_min_ =
          objective_function_values_candidates_ptr_[0];
      objective_function_values_min_index_ = 0;
      for (i = 1; i < pruned_set_size_; i++) {
        if (objective_function_values_candidates_ptr_[i] <
            objective_function_values_min_) {
          objective_function_values_min_ =
              objective_function_values_candidates_ptr_[i];
          objective_function_values_min_index_ = i;
        }
      }
      objective_function_values_(t) = objective_function_values_min_;
      change_points_[t] = pruned_set_[objective_function_values_min_index_];

      pruned_left_n_elem_ = 0;
      for (i = 0; i < pruned_set_size_; i++) {
        if (objective_function_values_candidates_ptr_[i] <=
            objective_function_values_min_ + beta_ - pruning_coefficient_) {
          pruned_set_[pruned_left_n_elem_] = pruned_set_[i];
          pruned_left_n_elem_++;
        }
      }
      pruned_set_size_ = pruned_left_n_elem_;
      pruned_set_[pruned_set_size_] = t;
      pruned_set_size_++;
    }
  } else if (family_ == "mean" && cost_adjustment_ == "BIC") {
    double two_norm;
    unsigned int i, pi;

    for (t = 2; t <= data_n_rows_; t++) {
      for (i = 0; i < pruned_set_size_; i++) {
        two_norm = (data_c_ptr_[t] - data_c_ptr_[pruned_set_[i]]) *
                   (data_c_ptr_[t] - data_c_ptr_[pruned_set_[i]]);
        for (pi = 1; pi < parameters_count_; pi++) {
          two_norm += (data_c_ptr_[t + data_c_n_rows_ * pi] -
                       data_c_ptr_[pruned_set_[i] + data_c_n_rows_ * pi]) *
                      (data_c_ptr_[t + data_c_n_rows_ * pi] -
                       data_c_ptr_[pruned_set_[i] + data_c_n_rows_ * pi]);
        }
        objective_function_values_candidates_ptr_[i] =
            objective_function_values_[pruned_set_[i]] +
            ((data_c_ptr_[t + data_c_n_rows_ * parameters_count_] -
              data_c_ptr_[pruned_set_[i] +
                          data_c_n_rows_ * parameters_count_]) -
             two_norm / (t - pruned_set_[i])) /
                2.0 +
            beta_;
      }

      objective_function_values_min_ =
          objective_function_values_candidates_ptr_[0];
      objective_function_values_min_index_ = 0;
      for (i = 1; i < pruned_set_size_; i++) {
        if (objective_function_values_candidates_ptr_[i] <
            objective_function_values_min_) {
          objective_function_values_min_ =
              objective_function_values_candidates_ptr_[i];
          objective_function_values_min_index_ = i;
        }
      }
      objective_function_values_(t) = objective_function_values_min_;
      change_points_[t] = pruned_set_[objective_function_values_min_index_];

      pruned_left_n_elem_ = 0;
      for (i = 0; i < pruned_set_size_; i++) {
        if (objective_function_values_candidates_ptr_[i] <=
            objective_function_values_min_ + beta_ - pruning_coefficient_) {
          pruned_set_[pruned_left_n_elem_] = pruned_set_[i];
          pruned_left_n_elem_++;
        }
      }
      pruned_set_size_ = pruned_left_n_elem_;
      pruned_set_[pruned_set_size_] = t;
      pruned_set_size_++;
    }
  } else if (family_ == "variance" && cost_adjustment_ == "MBIC" &&
             data_n_dims_ == 1) {
    unsigned int i;

    for (t = 2; t <= data_n_rows_; t++) {
      for (i = 0; i < pruned_set_size_; i++) {
        unsigned int const segment_length = t - pruned_set_[i];
        double det_value = data_c_ptr_[t] - data_c_ptr_[pruned_set_[i]];
        if (det_value <= 0) {
          det_value = 1e-11;
        }
        objective_function_values_candidates_ptr_[i] =
            objective_function_values_[pruned_set_[i]] +
            (log(det_value / segment_length) * segment_length / 2.0) +
            std::log(segment_length) / 2.0 + beta_;
      }

      objective_function_values_min_ =
          objective_function_values_candidates_ptr_[0];
      objective_function_values_min_index_ = 0;
      for (i = 1; i < pruned_set_size_; i++) {
        if (objective_function_values_candidates_ptr_[i] <
            objective_function_values_min_) {
          objective_function_values_min_ =
              objective_function_values_candidates_ptr_[i];
          objective_function_values_min_index_ = i;
        }
      }
      objective_function_values_(t) = objective_function_values_min_;
      change_points_[t] = pruned_set_[objective_function_values_min_index_];

      pruned_left_n_elem_ = 0;
      for (i = 0; i < pruned_set_size_; i++) {
        if (objective_function_values_candidates_ptr_[i] <=
            objective_function_values_min_ + beta_ - pruning_coefficient_) {
          pruned_set_[pruned_left_n_elem_] = pruned_set_[i];
          pruned_left_n_elem_++;
        }
      }
      pruned_set_size_ = pruned_left_n_elem_;
      pruned_set_[pruned_set_size_] = t;
      pruned_set_size_++;
    }
  } else if (family_ == "variance" && cost_adjustment_ == "MBIC" &&
             data_n_dims_ > 1) {
    unsigned int i;

    for (t = 2; t <= data_n_rows_; t++) {
      for (i = 0; i < pruned_set_size_; i++) {
        unsigned int const segment_length = t - pruned_set_[i];
        mat covar = zeros<mat>(data_n_dims_, data_n_dims_);
        for (unsigned int di = 0; di < data_n_dims_; di++) {
          for (unsigned int dj = 0; dj < data_n_dims_; dj++) {
            covar(di, dj) =
                data_c_ptr_[t + data_c_n_rows_ * (di * data_n_dims_ + dj)] -
                data_c_ptr_[pruned_set_[i] +
                            data_c_n_rows_ * (di * data_n_dims_ + dj)];
          }
        }
        double det_value = det(covar / segment_length);
        if (segment_length < data_n_dims_) {
          unsigned int approximate_segment_start;
          unsigned int approximate_segment_end;
          if (pruned_set_[i] >= data_n_dims_) {
            approximate_segment_start = pruned_set_[i] - data_n_dims_;
          } else {
            approximate_segment_start = 0;
          }
          if (t - 1 < data_n_rows_ - data_n_dims_) {
            approximate_segment_end = t - 1 + data_n_dims_;
          } else {
            approximate_segment_end = data_n_rows_ - 1;
          }
          mat covar_approx = zeros<mat>(data_n_dims_, data_n_dims_);
          for (unsigned int di = 0; di < data_n_dims_; di++) {
            for (unsigned int dj = 0; dj < data_n_dims_; dj++) {
              covar_approx(di, dj) =
                  data_c_ptr_[approximate_segment_end + 1 +
                              data_c_n_rows_ * (di * data_n_dims_ + dj)] -
                  data_c_ptr_[approximate_segment_start +
                              data_c_n_rows_ * (di * data_n_dims_ + dj)];
            }
          }
          det_value = det(covar_approx / (approximate_segment_end -
                                          approximate_segment_start + 1));
        }
        objective_function_values_candidates_ptr_[i] =
            objective_function_values_[pruned_set_[i]] +
            (log(det_value) * segment_length / 2.0) +
            std::log(segment_length) / 2.0 + beta_;
      }

      objective_function_values_min_ =
          objective_function_values_candidates_ptr_[0];
      objective_function_values_min_index_ = 0;
      for (i = 1; i < pruned_set_size_; i++) {
        if (objective_function_values_candidates_ptr_[i] <
            objective_function_values_min_) {
          objective_function_values_min_ =
              objective_function_values_candidates_ptr_[i];
          objective_function_values_min_index_ = i;
        }
      }
      objective_function_values_(t) = objective_function_values_min_;
      change_points_[t] = pruned_set_[objective_function_values_min_index_];

      pruned_left_n_elem_ = 0;
      for (i = 0; i < pruned_set_size_; i++) {
        if (objective_function_values_candidates_ptr_[i] <=
            objective_function_values_min_ + beta_ - pruning_coefficient_) {
          pruned_set_[pruned_left_n_elem_] = pruned_set_[i];
          pruned_left_n_elem_++;
        }
      }
      pruned_set_size_ = pruned_left_n_elem_;
      pruned_set_[pruned_set_size_] = t;
      pruned_set_size_++;
    }
  } else if (family_ == "meanvariance" && cost_adjustment_ == "MBIC" &&
             data_n_dims_ == 1) {
    unsigned int i;

    for (t = 2; t <= data_n_rows_; t++) {
      for (i = 0; i < pruned_set_size_; i++) {
        unsigned int const segment_length = t - pruned_set_[i];
        double det_value = data_c_ptr_[t + data_c_n_rows_] -
                           data_c_ptr_[pruned_set_[i] + data_c_n_rows_] -
                           (data_c_ptr_[t] - data_c_ptr_[pruned_set_[i]]) *
                               (data_c_ptr_[t] - data_c_ptr_[pruned_set_[i]]) /
                               segment_length;
        if (det_value <= 0) {
          det_value = 1e-11;
        }
        objective_function_values_candidates_ptr_[i] =
            objective_function_values_[pruned_set_[i]] +
            (log(det_value / segment_length) * segment_length / 2.0) +
            std::log(segment_length) / 2.0 + beta_;
      }

      objective_function_values_min_ =
          objective_function_values_candidates_ptr_[0];
      objective_function_values_min_index_ = 0;
      for (i = 1; i < pruned_set_size_; i++) {
        if (objective_function_values_candidates_ptr_[i] <
            objective_function_values_min_) {
          objective_function_values_min_ =
              objective_function_values_candidates_ptr_[i];
          objective_function_values_min_index_ = i;
        }
      }
      objective_function_values_(t) = objective_function_values_min_;
      change_points_[t] = pruned_set_[objective_function_values_min_index_];

      pruned_left_n_elem_ = 0;
      for (i = 0; i < pruned_set_size_; i++) {
        if (objective_function_values_candidates_ptr_[i] <=
            objective_function_values_min_ + beta_ - pruning_coefficient_) {
          pruned_set_[pruned_left_n_elem_] = pruned_set_[i];
          pruned_left_n_elem_++;
        }
      }
      pruned_set_size_ = pruned_left_n_elem_;
      pruned_set_[pruned_set_size_] = t;
      pruned_set_size_++;
    }
  } else if (family_ == "meanvariance" && cost_adjustment_ == "MBIC" &&
             data_n_dims_ > 1) {
    unsigned int i;

    for (t = 2; t <= data_n_rows_; t++) {
      for (i = 0; i < pruned_set_size_; i++) {
        unsigned int const segment_length = t - pruned_set_[i];
        mat covar = zeros<mat>(data_n_dims_, data_n_dims_);
        for (unsigned int di = 0; di < data_n_dims_; di++) {
          for (unsigned int dj = 0; dj < data_n_dims_; dj++) {
            unsigned int const crossprod_index =
                data_n_dims_ + di + data_n_dims_ * dj;
            covar(di, dj) =
                data_c_ptr_[t + data_c_n_rows_ * crossprod_index] -
                data_c_ptr_[pruned_set_[i] + data_c_n_rows_ * crossprod_index] -
                (data_c_ptr_[t + data_c_n_rows_ * di] -
                 data_c_ptr_[pruned_set_[i] + data_c_n_rows_ * di]) *
                    (data_c_ptr_[t + data_c_n_rows_ * dj] -
                     data_c_ptr_[pruned_set_[i] + data_c_n_rows_ * dj]) /
                    segment_length;
          }
        }
        double det_value = det(covar / segment_length);
        if (segment_length <= data_n_dims_) {
          unsigned int approximate_segment_start;
          unsigned int approximate_segment_end;
          if (pruned_set_[i] >= data_n_dims_) {
            approximate_segment_start = pruned_set_[i] - data_n_dims_;
          } else {
            approximate_segment_start = 0;
          }
          if (t - 1 < data_n_rows_ - data_n_dims_) {
            approximate_segment_end = t - 1 + data_n_dims_;
          } else {
            approximate_segment_end = data_n_rows_ - 1;
          }
          mat covar_approx = zeros<mat>(data_n_dims_, data_n_dims_);
          for (unsigned int di = 0; di < data_n_dims_; di++) {
            for (unsigned int dj = 0; dj < data_n_dims_; dj++) {
              unsigned int const crossprod_index =
                  data_n_dims_ + di + data_n_dims_ * dj;
              covar_approx(di, dj) =
                  data_c_ptr_[approximate_segment_end + 1 +
                              data_c_n_rows_ * crossprod_index] -
                  data_c_ptr_[approximate_segment_start +
                              data_c_n_rows_ * crossprod_index] -
                  (data_c_ptr_[approximate_segment_end + 1 +
                               data_c_n_rows_ * di] -
                   data_c_ptr_[approximate_segment_start +
                               data_c_n_rows_ * di]) *
                      (data_c_ptr_[approximate_segment_end + 1 +
                                   data_c_n_rows_ * dj] -
                       data_c_ptr_[approximate_segment_start +
                                   data_c_n_rows_ * dj]) /
                      (approximate_segment_end - approximate_segment_start + 1);
            }
          }
          det_value = det(covar_approx / (approximate_segment_end -
                                          approximate_segment_start + 1));
        }
        objective_function_values_candidates_ptr_[i] =
            objective_function_values_[pruned_set_[i]] +
            (log(det_value) * segment_length / 2.0) +
            std::log(segment_length) / 2.0 + beta_;
      }

      objective_function_values_min_ =
          objective_function_values_candidates_ptr_[0];
      objective_function_values_min_index_ = 0;
      for (i = 1; i < pruned_set_size_; i++) {
        if (objective_function_values_candidates_ptr_[i] <
            objective_function_values_min_) {
          objective_function_values_min_ =
              objective_function_values_candidates_ptr_[i];
          objective_function_values_min_index_ = i;
        }
      }
      objective_function_values_(t) = objective_function_values_min_;
      change_points_[t] = pruned_set_[objective_function_values_min_index_];

      pruned_left_n_elem_ = 0;
      for (i = 0; i < pruned_set_size_; i++) {
        if (objective_function_values_candidates_ptr_[i] <=
            objective_function_values_min_ + beta_ - pruning_coefficient_) {
          pruned_set_[pruned_left_n_elem_] = pruned_set_[i];
          pruned_left_n_elem_++;
        }
      }
      pruned_set_size_ = pruned_left_n_elem_;
      pruned_set_[pruned_set_size_] = t;
      pruned_set_size_++;
    }
  } else {
    CreateRProgress();
    UpdateRProgress();
    for (t = 2; t <= data_n_rows_; t++) {
      UpdateStep();
    }
  }
  return GetChangePointSet();
}

unordered_map<string,
              Fastcpd::FunctionSet> const Fastcpd::family_function_map_ = {
    {"arma", FunctionSet{&Fastcpd::GetGradientArma, &Fastcpd::GetHessianArma,
                         &Fastcpd::GetNllPeltArma, &Fastcpd::GetNllPeltArma,
                         &Fastcpd::GetNllSenArma}},
    {"binomial",
     FunctionSet{&Fastcpd::GetGradientBinomial, &Fastcpd::GetHessianBinomial,
                 &Fastcpd::GetNllPeltGlm, &Fastcpd::GetNllPeltGlm,
                 &Fastcpd::GetNllSenBinomial}},
    {"garch",
     FunctionSet{
         nullptr,  // No gradient
         nullptr,  // No hessian
         &Fastcpd::GetNllPeltGarch, &Fastcpd::GetNllPeltGarch,
         nullptr  // No nll_sen
     }},
    {"gaussian", FunctionSet{&Fastcpd::GetGradientLm, &Fastcpd::GetHessianLm,
                             &Fastcpd::GetNllPeltGlm, &Fastcpd::GetNllPeltGlm,
                             &Fastcpd::GetNllSenLm}},
    {"lasso", FunctionSet{&Fastcpd::GetGradientLm, &Fastcpd::GetHessianLm,
                          &Fastcpd::GetNllPeltLasso, &Fastcpd::GetNllPeltLasso,
                          &Fastcpd::GetNllSenLasso}},
    {"ma", FunctionSet{&Fastcpd::GetGradientMa, &Fastcpd::GetHessianMa,
                       &Fastcpd::GetNllPeltArma, &Fastcpd::GetNllPeltArma,
                       &Fastcpd::GetNllSenMa}},
    {"poisson",
     FunctionSet{&Fastcpd::GetGradientPoisson, &Fastcpd::GetHessianPoisson,
                 &Fastcpd::GetNllPeltGlm, &Fastcpd::GetNllPeltGlm,
                 &Fastcpd::GetNllSenPoisson}},
    {"mean",
     FunctionSet{
         nullptr,  // No gradient
         nullptr,  // No hessian
         &Fastcpd::GetNllPeltMean, &Fastcpd::GetNllPeltMeanValue,
         nullptr  // No nll_sen
     }},
    {"meanvariance",
     FunctionSet{
         nullptr,  // No gradient
         nullptr,  // No hessian
         &Fastcpd::GetNllPeltMeanvariance,
         &Fastcpd::GetNllPeltMeanvarianceValue,
         nullptr  // No nll_sen
     }},
    {"mgaussian",
     FunctionSet{
         nullptr,  // No gradient
         nullptr,  // No hessian
         &Fastcpd::GetNllPeltMgaussian, &Fastcpd::GetNllPeltMgaussian,
         nullptr  // No nll_sen
     }},
    {"variance",
     FunctionSet{
         nullptr,  // No gradient
         nullptr,  // No hessian
         &Fastcpd::GetNllPeltVariance, &Fastcpd::GetNllPeltVarianceValue,
         nullptr  // No nll_sen
     }}};

void Fastcpd::CreateRProgress() {
  if (r_progress_) {
#ifndef NO_RCPP
    rProgress_->tick(0);
#endif
  }
}

void Fastcpd::CreateSenParameters() {
  if (vanilla_percentage_ == 1) return;
  if (family_ == "binomial") {
    coefficients_sum_.col(0) = segment_coefficients_.row(0).t();
    coefficients_.col(0) = segment_coefficients_.row(0).t();
    double const prob =
        1 / (1 + exp(-dot(coefficients_.col(0),
                          data_.row(0).tail(parameters_count_))));
    hessian_.slice(0) = (data_.row(0).tail(parameters_count_).t() *
                         data_.row(0).tail(parameters_count_)) *
                        prob * (1 - prob);
  } else if (family_ == "poisson") {
    coefficients_.col(0) = segment_coefficients_.row(0).t();
    coefficients_sum_.col(0) = segment_coefficients_.row(0).t();
    hessian_.slice(0) =
        (data_.row(0).tail(parameters_count_).t() *
         data_.row(0).tail(parameters_count_)) *
        exp(dot(coefficients_.col(0), data_.row(0).tail(parameters_count_)));
  } else if (family_ == "lasso" || family_ == "gaussian") {
    coefficients_.col(0) = segment_coefficients_.row(0).t();
    coefficients_sum_.col(0) = segment_coefficients_.row(0).t();
    hessian_.slice(0) =
        data_.row(0).tail(parameters_count_).t() *
            data_.row(0).tail(parameters_count_) +
        epsilon_in_hessian_ * eye<mat>(parameters_count_, parameters_count_);
  } else if (family_ == "custom") {
    coefficients_.col(0) = segment_coefficients_.row(0).t();
    coefficients_sum_.col(0) = segment_coefficients_.row(0).t();
    hessian_.slice(0) = zeros<mat>(parameters_count_, parameters_count_);
  }
}

// TODO(doccstat): Use `segment_theta` as warm start.

void Fastcpd::CreateSegmentStatistics() {
  if (family_ == "mean" || family_ == "variance" || family_ == "meanvariance" ||
      (family_ == "custom" && vanilla_percentage_ == 1))
    return;
  for (int segment_index = 0; segment_index < segment_count_; ++segment_index) {
    GetCostResult(segment_indices_(segment_index),
                  segment_indices_(segment_index + 1) - 1, std::nullopt, true,
                  std::nullopt);

    // Initialize the estimated coefficients for each segment to be the
    // estimated coefficients in the segment.
    segment_coefficients_.row(segment_index) = result_coefficients_.t();
    if (family_ == "lasso" || family_ == "gaussian") {
      mat data_segment = data_.rows(segment_indices_(segment_index),
                                    segment_indices_(segment_index + 1) - 1);
      colvec segment_residual =
          data_segment.col(0) -
          data_segment.cols(1, data_segment.n_cols - 1) * result_coefficients_;
      double err_var = as_scalar(mean(square(segment_residual)));
      error_standard_deviation_(segment_index) = sqrt(err_var);
      active_coefficients_count_(segment_index) =
          accu(abs(result_coefficients_) > 0);
    }
  }
  // Mean of `error_standard_deviation_` only works if error sd is unchanged.
  lasso_penalty_base_ =
      mean(error_standard_deviation_) * sqrt(2 * std::log(parameters_count_));
  if (family_ == "lasso") {
    beta_ = beta_ * (1 + mean(active_coefficients_count_));
  }
}

double Fastcpd::GetCostAdjustmentValue(unsigned int const nrows) {
  double adjusted = 0;
  if (cost_adjustment_ == "MBIC" || cost_adjustment_ == "MDL") {
    adjusted = parameters_count_ * std::log(nrows) / 2.0;
  }
  if (cost_adjustment_ == "MDL") {
    adjusted *= std::log2(M_E);
  }
  return adjusted;
}

void Fastcpd::GetCostResult(unsigned int const segment_start,
                            unsigned int const segment_end,
                            std::optional<colvec> theta, bool const cv,
                            std::optional<colvec> start) {
  if (!theta.has_value()) {
    if (((family_ == "mean" || family_ == "variance" ||
          family_ == "meanvariance") &&
         t < data_c_n_rows_) ||
        (!(family_ == "mean" || family_ == "variance" ||
           family_ == "meanvariance") &&
         vanilla_percentage_ == 1)) {
      (this->*get_nll_pelt_value_)(segment_start, segment_end, cv, start);
    } else {
      (this->*get_nll_pelt_)(segment_start, segment_end, cv, start);
    }
  } else {
    result_coefficients_ = colvec();
    result_residuals_ = mat();
    colvec const theta_ = theta.value();
    result_value_ = (this->*get_nll_sen_)(segment_start, segment_end, theta_);
  }
  if (cost_adjustment_ == "MDL") {
    result_value_ = result_value_ * std::log2(M_E);
  }
  result_value_ += GetCostAdjustmentValue(segment_end - segment_start + 1);
}

std::tuple<arma::colvec, arma::colvec, arma::colvec, arma::mat, arma::mat>
Fastcpd::GetChangePointSet() {
  colvec cp_set = UpdateChangePointSet();

  if (cp_only_) {
    return std::make_tuple(change_points_, cp_set, colvec(), mat(), mat());
  }

  colvec cp_loc_ = zeros<colvec>(cp_set.n_elem + 2);
  if (cp_set.n_elem) {
    cp_loc_.rows(1, cp_loc_.n_elem - 2) = cp_set;
  }
  cp_loc_(cp_loc_.n_elem - 1) = data_n_rows_;
  colvec cp_loc = unique(std::move(cp_loc_));
  colvec cost_values = zeros<vec>(cp_loc.n_elem - 1);
  mat thetas = zeros<mat>(parameters_count_, cp_loc.n_elem - 1);
  mat residual;
  if (family_ == "mean" || family_ == "variance" || family_ == "meanvariance") {
    residual = zeros<mat>(data_n_rows_, data_n_dims_);
  } else if (family_ == "mgaussian") {
    residual = zeros<mat>(data_n_rows_, regression_response_count_);
  } else {
    residual = zeros<mat>(data_n_rows_, 1);
  }
  unsigned int residual_next_start = 0;

  for (unsigned int i = 0; i < cp_loc.n_elem - 1; i++) {
    GetCostResult(cp_loc(i), cp_loc(i + 1) - 1, std::nullopt, false,
                  std::nullopt);
    cost_values(i) = result_value_;
    if (family_ != "custom") {
      thetas.col(i) = result_coefficients_;
      residual.rows(residual_next_start,
                    residual_next_start + result_residuals_.n_rows - 1) =
          result_residuals_;
      residual_next_start += result_residuals_.n_rows;
    }
  }
  return std::make_tuple(change_points_, cp_set, cost_values, residual, thetas);
}

double Fastcpd::GetCostValue(int const tau, unsigned int const i) {
  if (t > vanilla_percentage_ * data_n_rows_) {
    return GetCostValueSen(tau, t - 1, i);
  } else {
    GetCostValuePelt(tau, t - 1, i);
    return result_value_;
  }
}

void Fastcpd::GetCostValuePelt(unsigned int const segment_start,
                               unsigned int const segment_end,
                               unsigned int const i) {
  if ((family_ == "binomial" || family_ == "poisson") &&
      (use_warm_start_ &&
       segment_end + 1 - segment_start >= 10 * parameters_count_)) {
    GetCostResult(segment_start, segment_end, std::nullopt, false,
                  segment_coefficients_
                      .row(index_max(find(segment_indices_ <= segment_end)))
                      .t()
                  // Or use `wrap(start.col(segment_start))` for warm start.
    );
    warm_start_.col(segment_start) = result_coefficients_;
  } else {
    GetCostResult(segment_start, segment_end, std::nullopt, false,
                  std::nullopt);
  }

  // If `vanilla_percentage_` is not 1, then we need to keep track of
  // thetas for later `fastcpd` steps.
  if (vanilla_percentage_ < 1 &&
      segment_end < vanilla_percentage_ * data_n_rows_) {
    coefficients_.col(i) = result_coefficients_;
    coefficients_sum_.col(i) += result_coefficients_;
  }
}

double Fastcpd::GetCostValueSen(unsigned int const segment_start,
                                unsigned int const segment_end,
                                unsigned int const i) {
  unsigned int const segment_length = segment_end - segment_start + 1;
  double cval = 0;
  UpdateSenParametersSteps(segment_start, segment_end, i);
  colvec theta = coefficients_sum_.col(i) / segment_length;
  if (family_ == "custom") {
    cval = (this->*get_nll_sen_)(segment_start, segment_end, theta);
  } else if ((family_ != "lasso" && segment_length >= parameters_count_) ||
             (family_ == "lasso" && segment_length >= 3)) {
    GetCostResult(segment_start, segment_end, theta, false, std::nullopt);
    cval = result_value_;
  }
  // else segment_length < parameters_count_ or for lasso segment_length < 3
  return cval;
}

#ifndef NO_RCPP
void Fastcpd::GetOptimizedCostResult(unsigned int const segment_start,
                                     unsigned int const segment_end) {
  mat const data_segment = data_.rows(segment_start, segment_end);
  if (parameters_count_ == 1) {
    Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
    Rcpp::Function optim = stats["optim"];
    Rcpp::List optim_result = optim(
        Rcpp::Named("par") = 0,
        Rcpp::Named("fn") = Rcpp::InternalFunction(
            +[](double theta, mat data_, Rcpp::Function cost) {
              return cost(Rcpp::Named("data") = data_,
                          Rcpp::Named("theta") = std::log(theta / (1 - theta)));
            }),
        Rcpp::Named("method") = "Brent", Rcpp::Named("lower") = 0,
        Rcpp::Named("upper") = 1, Rcpp::Named("data") = data_segment,
        Rcpp::Named("cost") = cost_function_.value());
    colvec par = Rcpp::as<colvec>(optim_result["par"]);
    double value = Rcpp::as<double>(optim_result["value"]);
    result_coefficients_ = log(par / (1 - par));
    result_residuals_ = mat();
    result_value_ = exp(value) / (1 + exp(value));
  } else {
    Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
    Rcpp::Function optim = stats["optim"];
    Rcpp::List optim_result = optim(
        Rcpp::Named("par") = zeros<vec>(parameters_count_),
        Rcpp::Named("fn") = cost_function_.value(),
        Rcpp::Named("method") = "L-BFGS-B", Rcpp::Named("data") = data_segment,
        Rcpp::Named("lower") = parameters_lower_bound_,
        Rcpp::Named("upper") = parameters_upper_bound_);
    result_coefficients_ = Rcpp::as<colvec>(optim_result["par"]);
    result_residuals_ = mat();
    result_value_ = Rcpp::as<double>(optim_result["value"]);
  }
}
#else
void Fastcpd::GetOptimizedCostResult(unsigned int const segment_start,
                                     unsigned int const segment_end) {
  // TODO(doccstat): Implement this function without Rcpp.
  // mat const data_segment = data_.rows(segment_start, segment_end);

  // optim::algo_settings_t settings;
  // settings.vals_bound = true;

  // colvec par = zeros<colvec>(parameters_count_);
  // settings.lower_bounds = parameters_lower_bound_;
  // settings.upper_bounds = parameters_upper_bound_;

  // auto objfn = [&](colvec const& x, colvec* /*grad_out*/,
  //                  void* /*opt_data*/) -> double {
  //   return cost_function_sen_(data_segment, x);
  // };
  // optim::lbfgs(par, objfn, nullptr, settings);
  // result_coefficients_ = std::move(par);
  // result_residuals_ = mat();
  // result_value_ = settings.opt_fn_value;
}
#endif

void Fastcpd::UpdateSenParametersStep(int const segment_start,
                                      int const segment_end, int const i) {
  mat hessian_i = hessian_.slice(i);
  colvec gradient;

  hessian_i +=
      (this->*get_hessian_)(segment_start, segment_end, coefficients_.col(i));
  gradient =
      (this->*get_gradient_)(segment_start, segment_end, coefficients_.col(i));

  // Add epsilon to the diagonal for PSD hessian_
  mat hessian_psd =
      hessian_i + epsilon_in_hessian_ *
                      eye<mat>(coefficients_.n_rows, coefficients_.n_rows);
  momentum_ = momentum_coef_ * momentum_ - solve(hessian_psd, gradient);
  double best_learning_rate = 1;
  colvec line_search_costs = zeros<colvec>(line_search_.n_elem);

  // Line search
  if (line_search_.n_elem > 1 || line_search_[0] != 1) {
    for (unsigned int line_search_index = 0;
         line_search_index < line_search_.n_elem; line_search_index++) {
      colvec theta_candidate =
          coefficients_.col(i) + line_search_[line_search_index] * momentum_;
      colvec theta_upper_bound =
          arma::min(std::move(theta_candidate), parameters_upper_bound_);
      colvec theta_projected =
          arma::max(std::move(theta_upper_bound), parameters_lower_bound_);
      line_search_costs[line_search_index] =
          (this->*get_nll_sen_)(segment_start, segment_end, theta_projected);
    }
  }
  best_learning_rate = line_search_[line_search_costs.index_min()];
  coefficients_.col(i) += best_learning_rate * momentum_;
  coefficients_.col(i) =
      arma::min(coefficients_.col(i), parameters_upper_bound_);
  coefficients_.col(i) =
      arma::max(coefficients_.col(i), parameters_lower_bound_);

  if (family_ == "lasso" || family_ == "gaussian") {
    // Update coefficients_ with L1 penalty
    double hessian_norm = norm(hessian_i, "fro");
    vec normd = abs(coefficients_.col(i));
    if (family_ == "lasso") {
      normd -= lasso_penalty_base_ / sqrt(segment_end - segment_start + 1) /
               hessian_norm;
    }
    coefficients_.col(i) = sign(coefficients_.col(i)) %
                           arma::max(normd, zeros<colvec>(normd.n_elem));
  }

  hessian_.slice(i) = std::move(hessian_i);
}

void Fastcpd::UpdateSenParametersSteps(int const segment_start,
                                       unsigned int const segment_end,
                                       int const i) {
  // This hack is to avoid the issue during momentum assigment of `solve`.
  colvec tmp = momentum_;
  unsigned int const multiple_epochs =
      multiple_epochs_function_(segment_end - segment_start + 1);
  unsigned int loop_start = segment_end, loop_end = segment_end;

  for (unsigned int epoch = 0; epoch <= multiple_epochs; epoch++) {
    for (loop_end = loop_start; loop_end <= segment_end; loop_end++) {
      UpdateSenParametersStep(segment_start, loop_end, i);
    }
    loop_start = segment_start;
  }

  coefficients_sum_.col(i) += coefficients_.col(i);
  momentum_ = tmp;
}

colvec Fastcpd::UpdateChangePointSet() {
  // Remove change points close to the boundaries.
  colvec cp_set = zeros<colvec>(data_n_rows_);
  int ncpts = 0;
  int last = data_n_rows_;
  while (last != 0) {
    cp_set[ncpts] = last;
    last = change_points_[last];
    ncpts += 1;
  }
  cp_set = sort(cp_set.rows(find(cp_set > 0)));
  cp_set = cp_set(find(cp_set > trim_ * data_n_rows_));
  cp_set = cp_set(find(cp_set < (1 - trim_) * data_n_rows_));
  colvec cp_set_ = zeros<vec>(cp_set.n_elem + 1);
  if (cp_set.n_elem) {
    cp_set_.rows(1, cp_set_.n_elem - 1) = std::move(cp_set);
  }
  cp_set = sort(unique(std::move(cp_set_)));

  // Remove change points close to each other.
  ucolvec cp_set_too_close = find(diff(cp_set) <= trim_ * data_n_rows_);
  if (cp_set_too_close.n_elem > 0) {
    int rest_element_count = cp_set.n_elem - cp_set_too_close.n_elem;
    colvec cp_set_rest_left = zeros<vec>(rest_element_count),
           cp_set_rest_right = zeros<vec>(rest_element_count);
    for (unsigned int i = 0, i_left = 0, i_right = 0; i < cp_set.n_elem; i++) {
      if (ucolvec left_find = find(cp_set_too_close == i);
          left_find.n_elem == 0) {
        cp_set_rest_left(i_left) = cp_set(i);
        i_left++;
      }
      if (ucolvec right_find = find(cp_set_too_close == i - 1);
          right_find.n_elem == 0) {
        cp_set_rest_right(i_right) = cp_set(i);
        i_right++;
      }
    }
    cp_set = floor((cp_set_rest_left + cp_set_rest_right) / 2);
  }
  return cp_set(find(cp_set > 0));
}

void Fastcpd::UpdateSenParameters() {
  if (vanilla_percentage_ == 1) return;
  int const segment_index = index_max(find(segment_indices_ <= t - 1));
  colvec cum_coef_add = segment_coefficients_.row(segment_index).t(),
         coef_add = segment_coefficients_.row(segment_index).t();
  mat hessian_new;
  if (family_ == "binomial") {
    rowvec const x = data_.row(t - 1).tail(parameters_count_);
    double const prob = 1 / (1 + exp(-dot(coef_add, x)));
    hessian_new = (x.t() * x) * prob * (1 - prob);
  } else if (family_ == "poisson") {
    rowvec const x = data_.row(t - 1).tail(parameters_count_);
    cum_coef_add = coef_add;
    hessian_new = (x.t() * x) * exp(dot(coef_add, x));
  } else if (family_ == "lasso" || family_ == "gaussian") {
    rowvec const x = data_.row(t - 1).tail(parameters_count_);
    hessian_new = x.t() * x + epsilon_in_hessian_ * eye<mat>(parameters_count_,
                                                             parameters_count_);
  } else if (family_ == "arma" || family_ == "ma") {
    hessian_new = (this->*get_hessian_)(0, t - 1, coef_add);
  } else if (family_ == "custom") {
    hessian_new = zeros<mat>(parameters_count_, parameters_count_);
  }
  memcpy(coefficients_.colptr(pruned_set_size_ - 1), coef_add.memptr(),
         sizeof(double) * parameters_count_);
  memcpy(coefficients_sum_.colptr(pruned_set_size_ - 1), cum_coef_add.memptr(),
         sizeof(double) * parameters_count_);
  memcpy(hessian_.slice(pruned_set_size_ - 1).memptr(), hessian_new.memptr(),
         sizeof(double) * parameters_count_ * parameters_count_);
}

void Fastcpd::UpdateRProgress() {
  if (r_progress_) {
#ifndef NO_RCPP
    rProgress_->tick();
#endif
  }
}

void Fastcpd::UpdateStep() {
  UpdateSenParameters();
  for (unsigned int i = 0; i < pruned_set_size_; i++) {
    if (i == pruned_set_size_ - 1 && vanilla_percentage_ != 1) {
      objective_function_values_candidates_ptr_[i] =
          objective_function_values_(pruned_set_(i)) + beta_;
    } else {
      objective_function_values_candidates_ptr_[i] =
          objective_function_values_(pruned_set_(i)) +
          GetCostValue(pruned_set_(i), i) + beta_;
    }
  }

  // The following code is the manual implementation of `index_min` function
  // in Armadillo.
  objective_function_values_min_ = objective_function_values_candidates_ptr_[0];
  objective_function_values_min_index_ = 0;
  for (unsigned int i = 1; i < pruned_set_size_; i++) {
    if (objective_function_values_candidates_ptr_[i] <
        objective_function_values_min_) {
      objective_function_values_min_ =
          objective_function_values_candidates_ptr_[i];
      objective_function_values_min_index_ = i;
    }
  }
  objective_function_values_(t) = objective_function_values_min_;
  change_points_[t] = pruned_set_[objective_function_values_min_index_];

  // The following code is the manual implementation of `find` function in
  // Armadillo.
  pruned_left_n_elem_ = 0;
  for (unsigned int i = 0; i < pruned_set_size_; i++) {
    if (objective_function_values_candidates_ptr_[i] <=
        objective_function_values_min_ + beta_ - pruning_coefficient_) {
      pruned_set_[pruned_left_n_elem_] = pruned_set_[i];
      if (vanilla_percentage_ != 1 && pruned_left_n_elem_ != i) {
        memcpy(coefficients_.colptr(pruned_left_n_elem_),
               coefficients_.colptr(i), sizeof(double) * parameters_count_);
        memcpy(coefficients_sum_.colptr(pruned_left_n_elem_),
               coefficients_sum_.colptr(i), sizeof(double) * parameters_count_);
        memcpy(hessian_.slice(pruned_left_n_elem_).memptr(),
               hessian_.slice(i).memptr(),
               sizeof(double) * parameters_count_ * parameters_count_);
      }
      pruned_left_[pruned_left_n_elem_] = i;
      pruned_left_n_elem_++;
    }
  }
  pruned_set_size_ = pruned_left_n_elem_;
  pruned_set_[pruned_set_size_] = t;
  pruned_set_size_++;
  UpdateRProgress();
}

colvec Fastcpd::GetGradientArma(unsigned int const segment_start,
                                unsigned int const segment_end,
                                colvec const& theta) {
  mat const data_segment = data_.rows(segment_start, segment_end);
  unsigned int const segment_length = segment_end - segment_start + 1;
  mat reversed_data = reverse(data_segment, 0);
  colvec reversed_theta = reverse(theta);
  if (segment_length < max(order_) + 1) {
    return ones(theta.n_elem);
  }
  colvec variance_term = zeros(segment_length);
  for (unsigned int i = max(order_); i < segment_length; i++) {
    variance_term(i) = data_segment(i, 0) -
                       dot(reversed_theta.rows(order_(1) + 1, sum(order_)),
                           data_segment.rows(i - order_(0), i - 1)) -
                       dot(reversed_theta.rows(1, order_(1)),
                           variance_term.rows(i - order_(1), i - 1));
  }
  colvec reversed_variance_term = reverse(variance_term);
  mat phi_coefficient = zeros(segment_length, order_(0)),
      psi_coefficient = zeros(segment_length, order_(1));
  for (unsigned int i = max(order_); i < segment_length; i++) {
    phi_coefficient.row(i) =
        -reversed_data
             .rows(segment_length - i, segment_length - i + order_(0) - 1)
             .t() -
        reversed_theta.rows(1, order_(1)).t() *
            phi_coefficient.rows(i - order_(1), i - 1);
  }
  for (unsigned int i = order_(1); i < segment_length; i++) {
    psi_coefficient.row(i) =
        -reversed_variance_term
             .rows(segment_length - i, segment_length - i + order_(1) - 1)
             .t() -
        reversed_theta.rows(1, order_(1)).t() *
            psi_coefficient.rows(i - order_(1), i - 1);
  }
  colvec gradient = zeros(sum(order_) + 1);
  gradient.rows(0, order_(0) - 1) =
      phi_coefficient.row(segment_length - 1).t() *
      variance_term(segment_length - 1) / theta(sum(order_));
  gradient.rows(order_(0), sum(order_) - 1) =
      psi_coefficient.row(segment_length - 1).t() *
      variance_term(segment_length - 1) / theta(sum(order_));
  gradient(sum(order_)) = 1.0 / 2.0 / theta(sum(order_)) -
                          pow(variance_term(segment_length - 1), 2) / 2.0 /
                              pow(theta(sum(order_)), 2);
  return gradient;
}

colvec Fastcpd::GetGradientBinomial(unsigned int const segment_start,
                                    unsigned int const segment_end,
                                    colvec const& theta) {
  mat const data_segment = data_.rows(segment_start, segment_end);
  unsigned int const segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double y = new_data(0);
  return -(y - 1 / (1 + exp(-as_scalar(x * theta)))) * x.t();
}

colvec Fastcpd::GetGradientCustom(unsigned int const segment_start,
                                  unsigned int const segment_end,
                                  colvec const& theta) {
  return cost_gradient_(data_.rows(segment_start, segment_end), theta);
}

colvec Fastcpd::GetGradientLm(unsigned int const segment_start,
                              unsigned int const segment_end,
                              colvec const& theta) {
  mat const data_segment = data_.rows(segment_start, segment_end);
  unsigned int const segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double y = new_data(0);
  return -(y - as_scalar(x * theta)) * x.t();
}

colvec Fastcpd::GetGradientMa(unsigned int const segment_start,
                              unsigned int const segment_end,
                              colvec const& theta) {
  mat const data_segment = data_.rows(segment_start, segment_end);
  unsigned int const segment_length = segment_end - segment_start + 1;
  unsigned int const q = order_(1);
  mat reversed_data = reverse(data_segment, 0);
  colvec reversed_theta = reverse(theta);
  if (segment_length < q + 1) {
    return ones(theta.n_elem);
  }
  colvec variance_term = zeros(segment_length);
  for (unsigned int i = q; i < segment_length; i++) {
    variance_term(i) =
        data_segment(i, 0) -
        dot(reversed_theta.rows(1, q), variance_term.rows(i - q, i - 1));
  }
  colvec reversed_variance_term = reverse(variance_term);
  mat psi_coefficient = zeros(segment_length, q);
  for (unsigned int i = q; i < segment_length; i++) {
    psi_coefficient.row(i) =
        -reversed_variance_term
             .rows(segment_length - i, segment_length - i + q - 1)
             .t() -
        reversed_theta.rows(1, q).t() * psi_coefficient.rows(i - q, i - 1);
  }
  colvec gradient = zeros(q + 1);
  gradient.rows(0, q - 1) = psi_coefficient.row(segment_length - 1).t() *
                            variance_term(segment_length - 1) / theta(q);
  gradient(q) =
      1.0 / 2.0 / theta(q) -
      pow(variance_term(segment_length - 1), 2) / 2.0 / pow(theta(q), 2);
  return gradient;
}

colvec Fastcpd::GetGradientPoisson(unsigned int const segment_start,
                                   unsigned int const segment_end,
                                   colvec const& theta) {
  mat const data_segment = data_.rows(segment_start, segment_end);
  unsigned int const segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double y = new_data(0);
  return -(y - exp(as_scalar(x * theta))) * x.t();
}

mat Fastcpd::GetHessianArma(unsigned int const segment_start,
                            unsigned int const segment_end,
                            colvec const& theta) {
  mat const data_segment = data_.rows(segment_start, segment_end);
  unsigned int const segment_length = segment_end - segment_start + 1;
  // TODO(doccstat): Maybe we can store all these computations
  mat reversed_data = reverse(data_segment, 0);
  colvec reversed_theta = reverse(theta);
  if (segment_length < max(order_) + 1) {
    return eye(theta.n_elem, theta.n_elem);
  }
  colvec variance_term = zeros(segment_length);
  for (unsigned int i = max(order_); i < segment_length; i++) {
    variance_term(i) = data_segment(i, 0) -
                       dot(reversed_theta.rows(order_(1) + 1, sum(order_)),
                           data_segment.rows(i - order_(0), i - 1)) -
                       dot(reversed_theta.rows(1, order_(1)),
                           variance_term.rows(i - order_(1), i - 1));
  }
  colvec reversed_variance_term = reverse(variance_term);
  mat phi_coefficient = zeros(segment_length, order_(0)),
      psi_coefficient = zeros(segment_length, order_(1));
  for (unsigned int i = max(order_); i < segment_length; i++) {
    phi_coefficient.row(i) =
        -reversed_data
             .rows(segment_length - i, segment_length - i + order_(0) - 1)
             .t() -
        reversed_theta.rows(1, order_(1)).t() *
            phi_coefficient.rows(i - order_(1), i - 1);
  }
  for (unsigned int i = order_(1); i < segment_length; i++) {
    psi_coefficient.row(i) =
        -reversed_variance_term
             .rows(segment_length - i, segment_length - i + order_(1) - 1)
             .t() -
        reversed_theta.rows(1, order_(1)).t() *
            psi_coefficient.rows(i - order_(1), i - 1);
  }
  mat reversed_coef_phi = reverse(phi_coefficient, 0),
      reversed_coef_psi = reverse(psi_coefficient, 0);
  cube phi_psi_coefficient = zeros(order_(1), order_(0), segment_length),
       psi_psi_coefficient = zeros(order_(1), order_(1), segment_length);
  for (unsigned int i = order_(1); i < segment_length; i++) {
    mat phi_psi_coefficient_part = zeros(order_(1), order_(0)),
        psi_psi_coefficient_part = zeros(order_(1), order_(1));
    for (unsigned int j = 1; j <= order_(1); j++) {
      phi_psi_coefficient_part +=
          phi_psi_coefficient.slice(i - j) * theta(order_(0) - 1 + j);
    }
    phi_psi_coefficient.slice(i) =
        -reversed_coef_phi.rows(segment_length - i,
                                segment_length - i + order_(1) - 1) -
        phi_psi_coefficient_part;
    for (unsigned int j = 1; j <= order_(1); j++) {
      psi_psi_coefficient_part +=
          psi_psi_coefficient.slice(i - j) * theta(order_(0) - 1 + j);
    }
    psi_psi_coefficient.slice(i) =
        -reversed_coef_psi.rows(segment_length - i,
                                segment_length - i + order_(1) - 1) -
        reversed_coef_psi
            .rows(segment_length - i, segment_length - i + order_(1) - 1)
            .t() -
        psi_psi_coefficient_part;
  }
  mat hessian = zeros(sum(order_) + 1, sum(order_) + 1);
  hessian.submat(0, 0, order_(0) - 1, order_(0) - 1) =
      phi_coefficient.row(segment_length - 1).t() *
      phi_coefficient.row(segment_length - 1) / theta(sum(order_));
  hessian.submat(0, order_(0), order_(0) - 1, sum(order_) - 1) =
      (phi_psi_coefficient.slice(segment_length - 1).t() *
           variance_term(segment_length - 1) +
       phi_coefficient.row(segment_length - 1).t() *
           psi_coefficient.row(segment_length - 1)) /
      theta(sum(order_));
  hessian.submat(order_(0), 0, sum(order_) - 1, order_(0) - 1) =
      hessian.submat(0, order_(0), order_(0) - 1, sum(order_) - 1).t();
  hessian.submat(0, sum(order_), order_(0) - 1, sum(order_)) =
      -phi_coefficient.row(segment_length - 1).t() *
      variance_term(segment_length - 1) / theta(sum(order_)) /
      theta(sum(order_));
  hessian.submat(sum(order_), 0, sum(order_), order_(0) - 1) =
      hessian.submat(0, sum(order_), order_(0) - 1, sum(order_)).t();
  hessian.submat(order_(0), order_(0), sum(order_) - 1, sum(order_) - 1) =
      (psi_coefficient.row(segment_length - 1).t() *
           psi_coefficient.row(segment_length - 1) +
       psi_psi_coefficient.slice(segment_length - 1) *
           variance_term(segment_length - 1)) /
      theta(sum(order_));
  hessian.submat(order_(0), sum(order_), sum(order_) - 1, sum(order_)) =
      -psi_coefficient.row(segment_length - 1).t() *
      variance_term(segment_length - 1) / theta(sum(order_)) /
      theta(sum(order_));
  hessian.submat(sum(order_), order_(0), sum(order_), sum(order_) - 1) =
      hessian.submat(order_(0), sum(order_), sum(order_) - 1, sum(order_)).t();
  hessian(sum(order_), sum(order_)) =
      pow(variance_term(segment_length - 1), 2) / pow(theta(sum(order_)), 3) -
      1.0 / 2.0 / pow(theta(sum(order_)), 2);
  return hessian;
}

mat Fastcpd::GetHessianBinomial(unsigned int const segment_start,
                                unsigned int const segment_end,
                                colvec const& theta) {
  mat const data_segment = data_.rows(segment_start, segment_end);
  unsigned int const segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double prob = 1 / (1 + exp(-dot(x, theta)));
  return (x.t() * x) * ((1 - prob) * prob);
}

mat Fastcpd::GetHessianCustom(unsigned int const segment_start,
                              unsigned int const segment_end,
                              colvec const& theta) {
  return cost_hessian_(data_.rows(segment_start, segment_end), theta);
}

mat Fastcpd::GetHessianLm(unsigned int const segment_start,
                          unsigned int const segment_end, colvec const& theta) {
  mat const data_segment = data_.rows(segment_start, segment_end);
  unsigned int const segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  return x.t() * x;
}

mat Fastcpd::GetHessianMa(unsigned int const segment_start,
                          unsigned int const segment_end, colvec const& theta) {
  mat const data_segment = data_.rows(segment_start, segment_end);
  unsigned int const segment_length = segment_end - segment_start + 1;
  unsigned int const q = order_(1);
  // TODO(doccstat): Maybe we can store all these computations
  mat reversed_data = reverse(data_segment, 0);
  colvec reversed_theta = reverse(theta);
  if (segment_length < q + 1) {
    return eye(theta.n_elem, theta.n_elem);
  }
  colvec variance_term = zeros(segment_length);
  for (unsigned int i = q; i < segment_length; i++) {
    variance_term(i) =
        data_segment(i, 0) -
        dot(reversed_theta.rows(1, q), variance_term.rows(i - q, i - 1));
  }
  colvec reversed_variance_term = reverse(variance_term);
  mat psi_coefficient = zeros(segment_length, q);
  for (unsigned int i = q; i < segment_length; i++) {
    psi_coefficient.row(i) =
        -reversed_variance_term
             .rows(segment_length - i, segment_length - i + q - 1)
             .t() -
        reversed_theta.rows(1, q).t() * psi_coefficient.rows(i - q, i - 1);
  }
  mat reversed_coef_psi = reverse(psi_coefficient, 0);
  cube psi_psi_coefficient = zeros(q, q, segment_length);
  for (unsigned int i = q; i < segment_length; i++) {
    mat psi_psi_coefficient_part = zeros(q, q);
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
  mat hessian = zeros(q + 1, q + 1);
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

mat Fastcpd::GetHessianPoisson(unsigned int const segment_start,
                               unsigned int const segment_end,
                               colvec const& theta) {
  mat const data_segment = data_.rows(segment_start, segment_end);
  unsigned int const segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double prob = exp(as_scalar(x * theta));
  // Prevent numerical issues if `prob` is too large.
  return (x.t() * x) * std::min(as_scalar(prob), 1e10);
}

#ifndef NO_RCPP
void Fastcpd::GetNllPeltArma(unsigned int const segment_start,
                             unsigned int const segment_end, bool const cv,
                             std::optional<colvec> const& start) {
  mat const data_segment = data_.rows(segment_start, segment_end);
  Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
  Rcpp::Function arima = stats["arima"];

  try {
    Rcpp::List out = arima(
        Rcpp::Named("x") = data_segment.col(0),
        Rcpp::Named("order") =
            Rcpp::NumericVector::create(order_(0), 0, order_(1)),
        Rcpp::Named("method") = "ML", Rcpp::Named("include.mean") = false);

    result_coefficients_ = zeros<colvec>(sum(order_) + 1);
    result_coefficients_.rows(0, sum(order_) - 1) =
        Rcpp::as<colvec>(out["coef"]);
    result_coefficients_(sum(order_)) = Rcpp::as<double>(out["sigma2"]);
    result_residuals_ = mat(Rcpp::as<colvec>(out["residuals"]));
    result_value_ = -Rcpp::as<double>(out["loglik"]);
  } catch (std::exception const& e) {
    // Handle the error - use reasonable defaults
    Rcpp::warning("ARMA model fitting failed: %s", e.what());

    // Set default coefficients (zeros)
    result_coefficients_ = zeros<colvec>(sum(order_) + 1);

    // Use a high penalty value to discourage this segment
    result_value_ = data_segment.n_rows * 10.0;

    // Set residuals Rcpp::as the original data (conservative approach)
    result_residuals_ = data_segment.col(0);

    // Set a relatively high variance estimate
    result_coefficients_(sum(order_)) = arma::var(data_segment.col(0));
  }
}
#else
void Fastcpd::GetNllPeltArma(unsigned int const segment_start,
                             unsigned int const segment_end, bool const cv,
                             std::optional<colvec> const& start) {
  // TODO(doccstat) : Implement this function.
}
#endif

void Fastcpd::GetNllPeltCustom(unsigned int const segment_start,
                               unsigned int const segment_end, bool const cv,
                               std::optional<colvec> const& start) {
  if (cost_gradient_ || cost_hessian_) {
    GetOptimizedCostResult(segment_start, segment_end);
  } else {
    result_coefficients_ = colvec();
    result_residuals_ = mat();
    result_value_ = cost_function_pelt_(data_.rows(segment_start, segment_end));
  }
}

#ifndef NO_RCPP
void Fastcpd::GetNllPeltGarch(unsigned int const segment_start,
                              unsigned int const segment_end, bool const cv,
                              std::optional<colvec> const& start) {
  colvec series = data_.rows(segment_start, segment_end).col(0);
  Rcpp::List out = garch(series, order_);
  result_coefficients_ = Rcpp::as<colvec>(out["coef"]);
  result_residuals_ = mat(Rcpp::as<colvec>(out["residuals"]));
  result_value_ = Rcpp::as<double>(out["n.likeli"]);
}
#else
void Fastcpd::GetNllPeltGarch(unsigned int const segment_start,
                              unsigned int const segment_end, bool const cv,
                              std::optional<colvec> const& start) {
  // TODO(doccstat) : Implement this function.
}
#endif

#ifndef NO_RCPP
void Fastcpd::GetNllPeltGlm(unsigned int const segment_start,
                            unsigned int const segment_end, bool const cv,
                            std::optional<colvec> const& start) {
  mat const data_segment = data_.rows(segment_start, segment_end);
  colvec y = data_segment.col(0);
  Rcpp::List out;
  if (!start.has_value()) {
    mat x = data_segment.cols(1, data_segment.n_cols - 1);
    out = fastglm(x, y, family_);
  } else {
    mat x = data_segment.cols(1, data_segment.n_cols - 1);
    out = fastglm(x, y, family_, start);
  }
  result_coefficients_ = Rcpp::as<colvec>(out["coefficients"]);
  result_residuals_ = mat(Rcpp::as<colvec>(out["residuals"]));
  result_value_ = Rcpp::as<double>(out["deviance"]) / 2;
}
#else
void Fastcpd::GetNllPeltGlm(unsigned int const segment_start,
                            unsigned int const segment_end, bool const cv,
                            std::optional<colvec> const& start) {
  // TODO(doccstat) : Implement this function.
}
#endif

#ifndef NO_RCPP
void Fastcpd::GetNllPeltLasso(unsigned int const segment_start,
                              unsigned int const segment_end, bool const cv,
                              std::optional<colvec> const& start) {
  if (segment_start == segment_end) {
    result_coefficients_ = zeros<colvec>(data_.n_cols - 1);
    result_residuals_ = zeros<mat>(1, 1);
    result_value_ = 0;
  } else if (cv) {
    mat const data_segment = data_.rows(segment_start, segment_end);
    Rcpp::Environment glmnet = Rcpp::Environment::namespace_env("glmnet"),
                      stats = Rcpp::Environment::namespace_env("stats");
    Rcpp::Function cv_glmnet = glmnet["cv.glmnet"],
                   predict_glmnet = glmnet["predict.glmnet"],
                   deviance = stats["deviance"];
    Rcpp::List out =
        cv_glmnet(data_segment.cols(1, data_segment.n_cols - 1),
                  data_segment.col(0), Rcpp::Named("family") = "gaussian");
    colvec index_vec = Rcpp::as<colvec>(out["index"]),
           values = Rcpp::as<colvec>(deviance(out["glmnet.fit"]));
    Rcpp::S4 out_coef = predict_glmnet(
        out["glmnet.fit"], Rcpp::Named("s") = out["lambda.1se"],
        Rcpp::Named("type") = "coefficients", Rcpp::Named("exact") = false);
    colvec glmnet_i = Rcpp::as<colvec>(out_coef.slot("i"));
    colvec glmnet_x = Rcpp::as<colvec>(out_coef.slot("x"));
    result_coefficients_ = zeros<colvec>(data_segment.n_cols - 1);
    for (unsigned int i = 1; i < glmnet_i.n_elem; i++) {
      result_coefficients_(glmnet_i(i) - 1) = glmnet_x(i);
    }
    result_residuals_ = mat();
    result_value_ = values(index_vec(1) - 1);
  } else {
    mat const data_segment = data_.rows(segment_start, segment_end);
    Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats"),
                      glmnet = Rcpp::Environment::namespace_env("glmnet");
    Rcpp::Function deviance = stats["deviance"], glmnet_ = glmnet["glmnet"],
                   predict_glmnet = glmnet["predict.glmnet"];
    Rcpp::List out =
        glmnet_(data_segment.cols(1, data_segment.n_cols - 1),
                data_segment.col(0), Rcpp::Named("family") = "gaussian",
                Rcpp::Named("lambda") = lasso_penalty_base_ /
                                        sqrt(segment_end - segment_start + 1));
    Rcpp::S4 out_par = out["beta"];
    colvec par_i = Rcpp::as<colvec>(out_par.slot("i"));
    colvec par_x = Rcpp::as<colvec>(out_par.slot("x"));
    result_coefficients_ = zeros<colvec>(data_segment.n_cols - 1);
    for (unsigned int i = 0; i < par_i.n_elem; i++) {
      result_coefficients_(par_i(i)) = par_x(i);
    }
    double value = Rcpp::as<double>(deviance(out));
    colvec fitted_values = Rcpp::as<colvec>(predict_glmnet(
        out, data_segment.cols(1, data_segment.n_cols - 1),
        Rcpp::Named("s") =
            lasso_penalty_base_ / sqrt(segment_end - segment_start + 1)));
    result_residuals_ = data_segment.col(0) - fitted_values;
    result_value_ = value / 2;
  }
}
#else
void Fastcpd::GetNllPeltLasso(unsigned int const segment_start,
                              unsigned int const segment_end, bool const cv,
                              std::optional<colvec> const& start) {
  // TODO(doccstat) : Implement this function.
}
#endif

void Fastcpd::GetNllPeltMean(unsigned int const segment_start,
                             unsigned int const segment_end, bool const cv,
                             std::optional<colvec> const& start) {
  mat data_segment = data_.rows(segment_start, segment_end);
  result_coefficients_ = mean(data_segment, 0).t();
  result_residuals_ = data_segment.each_row() - result_coefficients_.t();
  GetNllPeltMeanValue(segment_start, segment_end, cv, start);
}

void Fastcpd::GetNllPeltMeanValue(unsigned int const segment_start,
                                  unsigned int const segment_end, bool const cv,
                                  std::optional<colvec> const& start) {
  double two_norm = 0;
  for (unsigned int i = 0; i < parameters_count_; i++) {
    two_norm += (data_c_ptr_[segment_end + 1 + i * data_c_n_rows_] -
                 data_c_ptr_[segment_start + i * data_c_n_rows_]) *
                (data_c_ptr_[segment_end + 1 + i * data_c_n_rows_] -
                 data_c_ptr_[segment_start + i * data_c_n_rows_]);
  }
  unsigned int const segment_length = segment_end - segment_start + 1;
  result_value_ =
      ((data_c_ptr_[segment_end + 1 + parameters_count_ * data_c_n_rows_] -
        data_c_ptr_[segment_start + parameters_count_ * data_c_n_rows_]) -
       two_norm / segment_length) /
      2.0;
}

void Fastcpd::GetNllPeltMeanvariance(unsigned int const segment_start,
                                     unsigned int const segment_end,
                                     bool const cv,
                                     std::optional<colvec> const& start) {
  mat data_segment = data_.rows(segment_start, segment_end);
  mat cov_matrix = cov(data_segment);
  result_coefficients_ = zeros<colvec>(parameters_count_);
  result_coefficients_.rows(0, data_n_dims_ - 1) = mean(data_segment, 0).t();
  result_coefficients_.rows(data_n_dims_, parameters_count_ - 1) =
      cov_matrix.as_col();
  result_residuals_ = data_segment.each_row() - mean(data_segment, 0);
  result_residuals_.each_row() /= sqrt(cov_matrix.diag()).t();
  GetNllPeltMeanvarianceValue(segment_start, segment_end, cv, start);
}

void Fastcpd::GetNllPeltMeanvarianceValue(unsigned int const segment_start,
                                          unsigned int const segment_end,
                                          bool const cv,
                                          std::optional<colvec> const& start) {
  unsigned int const segment_length = segment_end - segment_start + 1;
  rowvec const data_diff =
      data_c_.row(segment_end + 1) - data_c_.row(segment_start);
  double det_value =
      det((reshape(data_diff.subvec(data_n_dims_, parameters_count_ - 1),
                   data_n_dims_, data_n_dims_) -
           (data_diff.subvec(0, data_n_dims_ - 1)).t() *
               (data_diff.subvec(0, data_n_dims_ - 1)) / segment_length) /
          segment_length);
  if (det_value <= 0) {
    det_value = 1e-10;
  }
  result_value_ = log(det_value) * segment_length / 2.0;
}

void Fastcpd::GetNllPeltMgaussian(unsigned int const segment_start,
                                  unsigned int const segment_end, bool const cv,
                                  std::optional<colvec> const& start) {
  mat const data_segment = data_.rows(segment_start, segment_end);
  mat x =
      data_segment.cols(regression_response_count_, data_segment.n_cols - 1);
  mat y = data_segment.cols(0, regression_response_count_ - 1);
  mat x_t_x;

  if (data_segment.n_rows <=
      data_segment.n_cols - regression_response_count_ + 1) {
    x_t_x = eye<mat>(data_segment.n_cols - regression_response_count_,
                     data_segment.n_cols - regression_response_count_);
  } else {
    x_t_x = x.t() * x;
  }

  mat coefficients = solve(x_t_x, x.t()) * y;
  result_coefficients_ = coefficients.as_col();
  result_residuals_ = y - x * coefficients;
  double value = regression_response_count_ * std::log(2.0 * M_PI) +
                 log_det_sympd(variance_estimate_);
  value *= data_segment.n_rows;
  value += trace(
      solve(variance_estimate_, result_residuals_.t() * result_residuals_));
  result_value_ = value / 2;
}

void Fastcpd::GetNllPeltVariance(unsigned int const segment_start,
                                 unsigned int const segment_end, bool const cv,
                                 std::optional<colvec> const& start) {
  mat data_segment = data_.rows(segment_start, segment_end);
  mat covar = cov(data_segment);
  result_coefficients_ = covar.as_col();
  result_residuals_ = data_segment.each_row() / sqrt(covar.diag()).t();
  GetNllPeltVarianceValue(segment_start, segment_end, cv, start);
}

void Fastcpd::GetNllPeltVarianceValue(unsigned int const segment_start,
                                      unsigned int const segment_end,
                                      bool const cv,
                                      std::optional<colvec> const& start) {
  unsigned int approximate_segment_start = segment_start,
               approximate_segment_end = segment_end;
  if (approximate_segment_end - approximate_segment_start + 1 < data_n_dims_) {
    if (segment_end < data_n_rows_ - data_n_dims_) {
      approximate_segment_end = segment_end + data_n_dims_;
    } else {
      approximate_segment_end = data_n_rows_ - 1;
    }
    approximate_segment_start = approximate_segment_end - data_n_dims_;
  }
  unsigned int const segment_length =
      approximate_segment_end - approximate_segment_start + 1;
  double const det_value =
      det(arma::reshape(data_c_.row(approximate_segment_end + 1) -
                            data_c_.row(approximate_segment_start),
                        data_n_dims_, data_n_dims_) /
          segment_length);
  result_value_ = log(det_value) * segment_length / 2.0;
}

double Fastcpd::GetNllSenArma(unsigned int const segment_start,
                              unsigned int const segment_end,
                              colvec const& theta) {
  mat data_segment = data_.rows(segment_start, segment_end);
  colvec reversed_theta = reverse(theta);
  if (data_segment.n_rows < max(order_) + 1) {
    return 0;
  }
  colvec variance_term = zeros(data_segment.n_rows);
  for (unsigned int i = max(order_); i < data_segment.n_rows; i++) {
    variance_term(i) = data_segment(i, 0) -
                       dot(reversed_theta.rows(order_(1) + 1, sum(order_)),
                           data_segment.rows(i - order_(0), i - 1)) -
                       dot(reversed_theta.rows(1, order_(1)),
                           variance_term.rows(i - order_(1), i - 1));
  }
  return (std::log(2.0 * M_PI) + std::log(theta(sum(order_)))) *
             (data_segment.n_rows - 2) / 2.0 +
         dot(variance_term, variance_term) / 2.0 / theta(sum(order_));
}

double Fastcpd::GetNllSenBinomial(unsigned int const segment_start,
                                  unsigned int const segment_end,
                                  colvec const& theta) {
  mat data_segment = data_.rows(segment_start, segment_end);
  colvec y = data_segment.col(0);
  // Calculate negative log likelihood in binomial family
  mat x = data_segment.cols(1, data_segment.n_cols - 1);
  colvec u = x * theta;
  return accu(-y % u + arma::log(1 + exp(u)));
}

double Fastcpd::GetNllSenCustom(unsigned int const segment_start,
                                unsigned int const segment_end,
                                colvec const& theta) {
  return cost_function_sen_(data_.rows(segment_start, segment_end), theta);
}

double Fastcpd::GetNllSenLasso(unsigned int const segment_start,
                               unsigned int const segment_end,
                               colvec const& theta) {
  mat data_segment = data_.rows(segment_start, segment_end);
  colvec y = data_segment.col(0);
  mat x = data_segment.cols(1, data_segment.n_cols - 1);
  return accu(square(y - x * theta)) / 2 +
         lasso_penalty_base_ / sqrt(segment_end - segment_start + 1) *
             accu(abs(theta));
}

double Fastcpd::GetNllSenLm(unsigned int const segment_start,
                            unsigned int const segment_end,
                            colvec const& theta) {
  mat data_segment = data_.rows(segment_start, segment_end);
  colvec y = data_segment.col(0);
  mat x = data_segment.cols(1, data_segment.n_cols - 1);
  return accu(square(y - x * theta)) / 2;
}

double Fastcpd::GetNllSenMa(unsigned int const segment_start,
                            unsigned int const segment_end,
                            colvec const& theta) {
  mat data_segment = data_.rows(segment_start, segment_end);
  unsigned int const q = order_(1);
  colvec reversed_theta = reverse(theta);
  if (data_segment.n_rows < q + 1) {
    return 0;
  }
  colvec variance_term = zeros(data_segment.n_rows);
  for (unsigned int i = q; i < data_segment.n_rows; i++) {
    variance_term(i) =
        data_segment(i, 0) -
        dot(reversed_theta.rows(1, q), variance_term.rows(i - q, i - 1));
  }
  return (std::log(2.0 * M_PI) + std::log(theta(q))) *
             (data_segment.n_rows - 2) / 2.0 +
         dot(variance_term, variance_term) / 2.0 / theta(q);
}

double Fastcpd::GetNllSenPoisson(unsigned int const segment_start,
                                 unsigned int const segment_end,
                                 colvec const& theta) {
  mat data_segment = data_.rows(segment_start, segment_end);
  colvec y = data_segment.col(0);
  mat x = data_segment.cols(1, data_segment.n_cols - 1);
  colvec u = x * theta;
  colvec y_factorial(y.n_elem);
  for (unsigned int i = 0; i < y.n_elem; i++) {
    double log_factorial = 0;
    for (int j = 1; j <= y(i); ++j) {
      log_factorial += std::log(j);
    }
    y_factorial(i) = log_factorial;
  }
  return accu(-y % u + exp(u) + y_factorial);
}

}  // namespace fastcpd::classes
