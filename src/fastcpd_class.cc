#include "fastcpd_class.h"

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
using ::arma::mat;
using ::arma::max;
using ::arma::mean;
using ::arma::min;
using ::arma::norm;
using ::arma::rowvec;
using ::arma::sign;
using ::arma::solve;
using ::arma::sort;
using ::arma::square;
using ::arma::ucolvec;
using ::arma::unique;
using ::arma::zeros;
using ::Rcpp::as;
using ::Rcpp::checkUserInterrupt;
using ::Rcpp::Function;
using ::Rcpp::List;
using ::Rcpp::Named;
using ::Rcpp::Nullable;
using ::Rcpp::wrap;
using ::std::make_unique;
using ::std::memcpy;
using ::std::string;
using ::std::string_view;
using ::std::unique_ptr;
using ::std::unordered_map;
using ::std::vector;

using ::arma::vec;
using ::Rcpp::Environment;
using ::Rcpp::InternalFunction;

using ::Rcpp::Rcout;

namespace fastcpd::classes {

Fastcpd::Fastcpd(const double beta, const Nullable<Function> cost,
                 const string cost_adjustment,
                 const Nullable<Function> cost_gradient,
                 const Nullable<Function> cost_hessian, const bool cp_only,
                 const unsigned int d, const mat data, const double epsilon,
                 const string family,
                 const Nullable<Function> multiple_epochs_function,
                 const colvec line_search, const colvec lower,
                 const double momentum_coef, const colvec order, const int p,
                 const unsigned int p_response, const double pruning_coef,
                 const string r_clock, const bool r_progress,
                 const int segment_count, const double trim, const colvec upper,
                 const double vanilla_percentage, const mat variance_estimate,
                 const bool warm_start)
    : active_coefficients_count_(colvec(segment_count)),
      beta_(beta),
      change_points_(zeros<colvec>(data.n_rows + 1)),
      coefficients_(mat(p, data.n_rows + 1)),
      coefficients_sum_(mat(p, data.n_rows + 1)),
      cost_adjustment_(cost_adjustment),
      cost_function_([&]() -> unique_ptr<Function> {
        if (family == "custom") {
          return make_unique<Function>(cost);
        }
        return nullptr;
      }()),
      cost_gradient_([&]() -> unique_ptr<Function> {
        if (family == "custom" && cost_gradient.isNotNull()) {
          return make_unique<Function>(cost_gradient);
        }
        return nullptr;
      }()),
      cost_hessian_([&]() -> unique_ptr<Function> {
        if (family == "custom" && cost_hessian.isNotNull()) {
          return make_unique<Function>(cost_hessian);
        }
        return nullptr;
      }()),
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
      data_n_dims_(d),
      data_n_cols_(data_.n_cols),
      data_n_rows_(data_.n_rows),
      epsilon_in_hessian_(epsilon),
      error_standard_deviation_(colvec(segment_count)),
      family_(family),
      hessian_(cube(p, p, data.n_rows + 1)),
      line_search_(line_search),
      momentum_(colvec(p)),
      momentum_coef_(momentum_coef),
      multiple_epochs_function_([&]() -> unique_ptr<Function> {
        if (multiple_epochs_function.isNotNull()) {
          return make_unique<Function>(multiple_epochs_function);
        }
        return nullptr;
      }()),
      objective_function_values_(colvec(data_n_rows_ + 1)),
      order_(order),
      parameters_count_(p),
      parameters_lower_bound_(lower),
      parameters_upper_bound_(upper),
      pruned_left_(ucolvec(data_n_rows_ + 1)),
      pruned_set_(zeros<ucolvec>(data_n_rows_ + 1)),
      pruning_coefficient_(pruning_coef),
      r_clock_(r_clock),
      r_progress_(r_progress),
      regression_response_count_(p_response),
      rProgress_(make_unique<RProgress::RProgress>(kRProgress, data_n_rows_)),
      segment_coefficients_(mat(segment_count, p)),
      segment_count_(segment_count),
      segment_indices_(round(linspace(0, data_n_rows_, segment_count + 1))),
      trim_(trim),
      use_warm_start_(warm_start),
      vanilla_percentage_(vanilla_percentage),
      variance_estimate_(variance_estimate),
      warm_start_(zeros<mat>(p, data_n_rows_)) {
  // Handle the special 'arma' family_ with order_ condition
  if (family_ == "arma") {
    if (order_(0) > 0) {
      get_gradient_ = &Fastcpd::GetGradientArma;
      get_hessian_ = &Fastcpd::GetHessianArma;
      get_nll_sen_ = &Fastcpd::GetNllSenArma;
      get_nll_pelt_ = &Fastcpd::GetNllPeltArma;
    } else {  // order_(0) == 0
      get_gradient_ = &Fastcpd::GetGradientMa;
      get_hessian_ = &Fastcpd::GetHessianMa;
      get_nll_sen_ = &Fastcpd::GetNllSenMa;
      get_nll_pelt_ = &Fastcpd::GetNllPeltArma;
    }
  } else {
    auto it = family_function_map_.find(family_);
    if (it != family_function_map_.end()) {
      const FunctionSet &func_set = it->second;
      get_gradient_ = func_set.gradient;
      get_hessian_ = func_set.hessian;
      get_nll_sen_ = func_set.nll_sen;
      get_nll_pelt_ = func_set.nll_pelt;
    } else {
      get_gradient_ = &Fastcpd::GetGradientCustom;
      get_hessian_ = &Fastcpd::GetHessianCustom;
      get_nll_sen_ = &Fastcpd::GetNllSenCustom;
      get_nll_pelt_ = &Fastcpd::GetNllPeltCustom;
    }
  }

  // TODO(doccstat): Store environment functions from R.
}

List Fastcpd::Run() {
  pruned_set_(1) = 1;
  objective_function_values_.fill(arma::datum::inf);
  objective_function_values_(0) = -beta_;
  objective_function_values_(1) = GetCostValue(0, 0, 1);

  // TODO(doccstat): Investigate if the following branches can be merged into
  // `fastcpd_class_nll.cc`.
  if (family_ == "mean" && cost_adjustment_ == "MBIC") {
    double *obj = (double *)calloc(data_n_rows_ + 1, sizeof(double));
    double two_norm;
    unsigned int i, pi;

    for (unsigned int t = 2; t <= data_n_rows_; t++) {
      for (i = 0; i < pruned_set_size_; i++) {
        two_norm = (data_c_ptr_[t] - data_c_ptr_[pruned_set_[i]]) *
                   (data_c_ptr_[t] - data_c_ptr_[pruned_set_[i]]);
        for (pi = 1; pi < parameters_count_; pi++) {
          two_norm += (data_c_ptr_[t + data_c_n_rows_ * pi] -
                       data_c_ptr_[pruned_set_[i] + data_c_n_rows_ * pi]) *
                      (data_c_ptr_[t + data_c_n_rows_ * pi] -
                       data_c_ptr_[pruned_set_[i] + data_c_n_rows_ * pi]);
        }
        obj[i] = objective_function_values_[pruned_set_[i]] +
                 ((data_c_ptr_[t + data_c_n_rows_ * parameters_count_] -
                   data_c_ptr_[pruned_set_[i] +
                               data_c_n_rows_ * parameters_count_]) -
                  two_norm / (t - pruned_set_[i])) /
                     2.0 +
                 std::log(t - pruned_set_[i]) / 2.0 + beta_;
      }

      objective_function_values_min_ = obj[0];
      objective_function_values_min_index_ = 0;
      for (i = 1; i < pruned_set_size_; i++) {
        if (obj[i] < objective_function_values_min_) {
          objective_function_values_min_ = obj[i];
          objective_function_values_min_index_ = i;
        }
      }
      objective_function_values_(t) = objective_function_values_min_;
      change_points_[t] = pruned_set_[objective_function_values_min_index_];

      pruned_left_n_elem_ = 0;
      for (i = 0; i < pruned_set_size_; i++) {
        if (obj[i] <=
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
    double *obj = (double *)calloc(data_n_rows_ + 1, sizeof(double));
    double two_norm;
    unsigned int i, pi;

    for (unsigned int t = 2; t <= data_n_rows_; t++) {
      for (i = 0; i < pruned_set_size_; i++) {
        two_norm = (data_c_ptr_[t] - data_c_ptr_[pruned_set_[i]]) *
                   (data_c_ptr_[t] - data_c_ptr_[pruned_set_[i]]);
        for (pi = 1; pi < parameters_count_; pi++) {
          two_norm += (data_c_ptr_[t + data_c_n_rows_ * pi] -
                       data_c_ptr_[pruned_set_[i] + data_c_n_rows_ * pi]) *
                      (data_c_ptr_[t + data_c_n_rows_ * pi] -
                       data_c_ptr_[pruned_set_[i] + data_c_n_rows_ * pi]);
        }
        obj[i] = objective_function_values_[pruned_set_[i]] +
                 ((data_c_ptr_[t + data_c_n_rows_ * parameters_count_] -
                   data_c_ptr_[pruned_set_[i] +
                               data_c_n_rows_ * parameters_count_]) -
                  two_norm / (t - pruned_set_[i])) /
                     2.0 +
                 beta_;
      }

      objective_function_values_min_ = obj[0];
      objective_function_values_min_index_ = 0;
      for (i = 1; i < pruned_set_size_; i++) {
        if (obj[i] < objective_function_values_min_) {
          objective_function_values_min_ = obj[i];
          objective_function_values_min_index_ = i;
        }
      }
      objective_function_values_(t) = objective_function_values_min_;
      change_points_[t] = pruned_set_[objective_function_values_min_index_];

      pruned_left_n_elem_ = 0;
      for (i = 0; i < pruned_set_size_; i++) {
        if (obj[i] <=
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
    double *obj = (double *)calloc(data_n_rows_ + 1, sizeof(double));
    unsigned int i;

    for (unsigned int t = 2; t <= data_n_rows_; t++) {
      for (i = 0; i < pruned_set_size_; i++) {
        const unsigned int segment_length = t - pruned_set_[i];
        double det_value = data_c_ptr_[t] - data_c_ptr_[pruned_set_[i]];
        if (det_value <= 0) {
          det_value = 1e-11;
        }
        obj[i] = objective_function_values_[pruned_set_[i]] +
                 (log(det_value / segment_length) * segment_length / 2.0) +
                 std::log(segment_length) / 2.0 + beta_;
      }

      objective_function_values_min_ = obj[0];
      objective_function_values_min_index_ = 0;
      for (i = 1; i < pruned_set_size_; i++) {
        if (obj[i] < objective_function_values_min_) {
          objective_function_values_min_ = obj[i];
          objective_function_values_min_index_ = i;
        }
      }
      objective_function_values_(t) = objective_function_values_min_;
      change_points_[t] = pruned_set_[objective_function_values_min_index_];

      pruned_left_n_elem_ = 0;
      for (i = 0; i < pruned_set_size_; i++) {
        if (obj[i] <=
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
    double *obj = (double *)calloc(data_n_rows_ + 1, sizeof(double));
    unsigned int i;

    for (unsigned int t = 2; t <= data_n_rows_; t++) {
      for (i = 0; i < pruned_set_size_; i++) {
        const unsigned int segment_length = t - pruned_set_[i];
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
        obj[i] = objective_function_values_[pruned_set_[i]] +
                 (log(det_value) * segment_length / 2.0) +
                 std::log(segment_length) / 2.0 + beta_;
      }

      objective_function_values_min_ = obj[0];
      objective_function_values_min_index_ = 0;
      for (i = 1; i < pruned_set_size_; i++) {
        if (obj[i] < objective_function_values_min_) {
          objective_function_values_min_ = obj[i];
          objective_function_values_min_index_ = i;
        }
      }
      objective_function_values_(t) = objective_function_values_min_;
      change_points_[t] = pruned_set_[objective_function_values_min_index_];

      pruned_left_n_elem_ = 0;
      for (i = 0; i < pruned_set_size_; i++) {
        if (obj[i] <=
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
    double *obj = (double *)calloc(data_n_rows_ + 1, sizeof(double));
    unsigned int i;

    for (unsigned int t = 2; t <= data_n_rows_; t++) {
      for (i = 0; i < pruned_set_size_; i++) {
        const unsigned int segment_length = t - pruned_set_[i];
        double det_value = data_c_ptr_[t + data_c_n_rows_] -
                           data_c_ptr_[pruned_set_[i] + data_c_n_rows_] -
                           (data_c_ptr_[t] - data_c_ptr_[pruned_set_[i]]) *
                               (data_c_ptr_[t] - data_c_ptr_[pruned_set_[i]]) /
                               segment_length;
        if (det_value <= 0) {
          det_value = 1e-11;
        }
        obj[i] = objective_function_values_[pruned_set_[i]] +
                 (log(det_value / segment_length) * segment_length / 2.0) +
                 std::log(segment_length) / 2.0 + beta_;
      }

      objective_function_values_min_ = obj[0];
      objective_function_values_min_index_ = 0;
      for (i = 1; i < pruned_set_size_; i++) {
        if (obj[i] < objective_function_values_min_) {
          objective_function_values_min_ = obj[i];
          objective_function_values_min_index_ = i;
        }
      }
      objective_function_values_(t) = objective_function_values_min_;
      change_points_[t] = pruned_set_[objective_function_values_min_index_];

      pruned_left_n_elem_ = 0;
      for (i = 0; i < pruned_set_size_; i++) {
        if (obj[i] <=
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
    double *obj = (double *)calloc(data_n_rows_ + 1, sizeof(double));
    unsigned int i;

    for (unsigned int t = 2; t <= data_n_rows_; t++) {
      for (i = 0; i < pruned_set_size_; i++) {
        const unsigned int segment_length = t - pruned_set_[i];
        mat covar = zeros<mat>(data_n_dims_, data_n_dims_);
        for (unsigned int di = 0; di < data_n_dims_; di++) {
          for (unsigned int dj = 0; dj < data_n_dims_; dj++) {
            const unsigned int crossprod_index =
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
              const unsigned int crossprod_index =
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
        obj[i] = objective_function_values_[pruned_set_[i]] +
                 (log(det_value) * segment_length / 2.0) +
                 std::log(segment_length) / 2.0 + beta_;
      }

      objective_function_values_min_ = obj[0];
      objective_function_values_min_index_ = 0;
      for (i = 1; i < pruned_set_size_; i++) {
        if (obj[i] < objective_function_values_min_) {
          objective_function_values_min_ = obj[i];
          objective_function_values_min_index_ = i;
        }
      }
      objective_function_values_(t) = objective_function_values_min_;
      change_points_[t] = pruned_set_[objective_function_values_min_index_];

      pruned_left_n_elem_ = 0;
      for (i = 0; i < pruned_set_size_; i++) {
        if (obj[i] <=
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
    CreateSegmentStatistics();
    CreateSenParameters();
    CreateRProgress();
    UpdateRProgress();
    for (unsigned int t = 2; t <= data_n_rows_; t++) {
      UpdateStep(t);
    }
  }
  return GetChangePointSet();
}

const unordered_map<string, Fastcpd::FunctionSet>
    Fastcpd::family_function_map_ = {
        {"binomial",
         FunctionSet{&Fastcpd::GetGradientBinomial,
                     &Fastcpd::GetHessianBinomial, &Fastcpd::GetNllSenBinomial,
                     &Fastcpd::GetNllPeltGlm}},
        {"garch", FunctionSet{nullptr,  // No gradient
                              nullptr,  // No hessian
                              nullptr,  // No nll_sen
                              &Fastcpd::GetNllPeltGarch}},
        {"gaussian",
         FunctionSet{&Fastcpd::GetGradientLm, &Fastcpd::GetHessianLm,
                     &Fastcpd::GetNllSenLm, &Fastcpd::GetNllPeltGlm}},
        {"lasso",
         FunctionSet{&Fastcpd::GetGradientLm, &Fastcpd::GetHessianLm,
                     &Fastcpd::GetNllSenLasso, &Fastcpd::GetNllPeltLasso}},
        {"poisson",
         FunctionSet{&Fastcpd::GetGradientPoisson, &Fastcpd::GetHessianPoisson,
                     &Fastcpd::GetNllSenPoisson, &Fastcpd::GetNllPeltGlm}},
        {"mean", FunctionSet{nullptr,  // No gradient
                             nullptr,  // No hessian
                             nullptr,  // No nll_sen
                             &Fastcpd::GetNllPeltMean}},
        {"meanvariance", FunctionSet{nullptr,  // No gradient
                                     nullptr,  // No hessian
                                     nullptr,  // No nll_sen
                                     &Fastcpd::GetNllPeltMeanVariance}},
        {"mgaussian", FunctionSet{nullptr,  // No gradient
                                  nullptr,  // No hessian
                                  nullptr,  // No nll_sen
                                  &Fastcpd::GetNllPeltMgaussian}},
        {"variance", FunctionSet{nullptr,  // No gradient
                                 nullptr,  // No hessian
                                 nullptr,  // No nll_sen
                                 &Fastcpd::GetNllPeltVariance}}};

void Fastcpd::CreateRClock(const std::string name) {
  if (!r_clock_.empty()) {
    rClock_.stop(name);
  }
}

void Fastcpd::CreateRProgress() {
  if (r_progress_) {
    rProgress_->tick(0);
  }
}

void Fastcpd::CreateSenParameters() {
  if (vanilla_percentage_ == 1) return;
  if (family_ == "binomial") {
    coefficients_sum_.col(0) = segment_coefficients_.row(0).t();
    coefficients_.col(0) = segment_coefficients_.row(0).t();
    const double prob =
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
  if (family_ == "custom" && vanilla_percentage_ == 1) return;
  for (int segment_index = 0; segment_index < segment_count_; ++segment_index) {
    GetCostResult(segment_indices_(segment_index),
                  segment_indices_(segment_index + 1) - 1, R_NilValue, true,
                  R_NilValue);
    colvec segment_theta = result_coefficients_.as_col();

    // Initialize the estimated coefficients for each segment to be the
    // estimated coefficients in the segment.
    segment_coefficients_.row(segment_index) = segment_theta.t();
    if (family_ == "lasso" || family_ == "gaussian") {
      mat data_segment = data_.rows(segment_indices_(segment_index),
                                    segment_indices_(segment_index + 1) - 1);
      colvec segment_residual =
          data_segment.col(0) -
          data_segment.cols(1, data_segment.n_cols - 1) * segment_theta;
      double err_var = as_scalar(mean(square(segment_residual)));
      error_standard_deviation_(segment_index) = sqrt(err_var);
      active_coefficients_count_(segment_index) = accu(abs(segment_theta) > 0);
    }
  }
  // Mean of `error_standard_deviation_` only works if error sd is unchanged.
  lasso_penalty_base_ =
      mean(error_standard_deviation_) * sqrt(2 * std::log(parameters_count_));
  if (family_ == "lasso") {
    beta_ = beta_ * (1 + mean(active_coefficients_count_));
  }
}

double Fastcpd::GetCostAdjustmentValue(const unsigned nrows) {
  double adjusted = 0;
  if (cost_adjustment_ == "MBIC" || cost_adjustment_ == "MDL") {
    adjusted = parameters_count_ * std::log(nrows) / 2.0;
  }
  if (cost_adjustment_ == "MDL") {
    adjusted *= std::log2(M_E);
  }
  return adjusted;
}

void Fastcpd::GetCostResult(const unsigned int segment_start,
                            const unsigned int segment_end,
                            Nullable<colvec> theta, const bool cv,
                            Nullable<colvec> start) {
  if (theta.isNull()) {
    (this->*get_nll_pelt_)(segment_start, segment_end, cv, start);
  } else {
    result_coefficients_ = mat();
    result_residuals_ = mat();
    result_value_ =
        (this->*get_nll_sen_)(segment_start, segment_end, as<colvec>(theta));
  }
  result_value_ =
      UpdateCostValue(result_value_, segment_end - segment_start + 1);
}

List Fastcpd::GetChangePointSet() {
  CreateRClock(r_clock_);
  colvec cp_set = UpdateChangePointSet();

  if (cp_only_) {
    return List::create(
        Named("raw_cp_set") = change_points_, Named("cp_set") = cp_set,
        Named("cost_values") = R_NilValue, Named("residual") = R_NilValue,
        Named("thetas") = R_NilValue);
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
    GetCostResult(cp_loc(i), cp_loc(i + 1) - 1, R_NilValue, false, R_NilValue);
    cost_values(i) = result_value_;

    // Parameters are not involved for PELT.
    if (vanilla_percentage_ < 1 || family_ == "garch") {
      thetas.col(i) = result_coefficients_.as_col();
    }

    // Residual is only calculated for built-in families.
    if (family_ != "custom" && family_ != "mean" && family_ != "variance" &&
        family_ != "meanvariance") {
      mat cost_optim_residual = result_residuals_;
      residual.rows(residual_next_start,
                    residual_next_start + cost_optim_residual.n_rows - 1) =
          cost_optim_residual;
      residual_next_start += cost_optim_residual.n_rows;
    }
  }
  return List::create(Named("raw_cp_set") = change_points_,
                      Named("cp_set") = cp_set,
                      Named("cost_values") = cost_values,
                      Named("residual") = residual, Named("thetas") = thetas);
}

double Fastcpd::GetCostValue(const int tau, const unsigned int i, const int t) {
  if (t > vanilla_percentage_ * data_n_rows_) {
    return GetCostValueSen(tau, t - 1, i);
  } else {
    return GetCostValuePelt(tau, t - 1, i);
  }
}

double Fastcpd::GetCostValuePelt(const unsigned int segment_start,
                                 const unsigned int segment_end,
                                 const unsigned int i) {
  double cval = 0;
  if ((family_ == "binomial" || family_ == "poisson") &&
      (use_warm_start_ &&
       segment_end + 1 - segment_start >= 10 * parameters_count_)) {
    GetCostResult(
        segment_start, segment_end, R_NilValue, false,
        wrap(segment_coefficients_
                 .row(index_max(find(segment_indices_ <= segment_end)))
                 .t())
        // Or use `wrap(start.col(segment_start))` for warm start.
    );
    warm_start_.col(segment_start) = result_coefficients_.as_col();
  } else {
    GetCostResult(segment_start, segment_end, R_NilValue, false, R_NilValue);
  }
  cval = result_value_;

  // If `vanilla_percentage_` is not 1, then we need to keep track of
  // thetas for later `fastcpd` steps.
  if (vanilla_percentage_ < 1 &&
      segment_end < vanilla_percentage_ * data_n_rows_) {
    coefficients_.col(i) = result_coefficients_.as_col();
    coefficients_sum_.col(i) += result_coefficients_.as_col();
  }
  return cval;
}

double Fastcpd::GetCostValueSen(const unsigned int segment_start,
                                const unsigned int segment_end,
                                const unsigned int i) {
  const unsigned int segment_length = segment_end - segment_start + 1;
  double cval = 0;
  UpdateSenParametersSteps(segment_start, segment_end, i);
  colvec theta = coefficients_sum_.col(i) / segment_length;
  if (family_ == "custom") {
    cval = (this->*get_nll_sen_)(segment_start, segment_end, theta);
  } else if ((family_ != "lasso" && segment_length >= parameters_count_) ||
             (family_ == "lasso" && segment_length >= 3)) {
    GetCostResult(segment_start, segment_end, wrap(theta), false, R_NilValue);
    cval = result_value_;
  }
  // else segment_length < parameters_count_ or for lasso segment_length < 3
  return cval;
}

colvec Fastcpd::GetObjectiveFunctionValues(unsigned int t) {
  colvec cval = zeros<vec>(pruned_set_size_);
  UpdateRClockTick("r_t_set_for_loop");
  unsigned int loop_end = pruned_set_size_ - (vanilla_percentage_ != 1);
  for (unsigned int i = 0; i < loop_end; i++) {
    cval(i) = GetCostValue(pruned_set_(i), i, t);
  }
  UpdateRClockTock("r_t_set_for_loop");
  colvec obj = cval +
               objective_function_values_.rows(
                   pruned_set_.rows(0, pruned_set_size_ - 1)) +
               beta_;
  return obj;
}

void Fastcpd::GetOptimizedCostResult(const unsigned int segment_start,
                                     const unsigned int segment_end) {
  const mat data_segment = data_.rows(segment_start, segment_end);
  if (parameters_count_ == 1) {
    Environment stats = Environment::namespace_env("stats");
    Function optim = stats["optim"];
    List optim_result =
        optim(Named("par") = 0,
              Named("fn") =
                  InternalFunction(+[](double theta, mat data_, Function cost) {
                    return cost(Named("data") = data_,
                                Named("theta") = std::log(theta / (1 - theta)));
                  }),
              Named("method") = "Brent", Named("lower") = 0, Named("upper") = 1,
              Named("data") = data_segment, Named("cost") = *cost_function_);
    colvec par = as<colvec>(optim_result["par"]);
    double value = as<double>(optim_result["value"]);
    result_coefficients_ = mat(log(par / (1 - par)));
    result_residuals_ = mat();
    result_value_ = exp(value) / (1 + exp(value));
  } else {
    Environment stats = Environment::namespace_env("stats");
    Function optim = stats["optim"];
    List optim_result = optim(
        Named("par") = zeros<vec>(parameters_count_),
        Named("fn") = *cost_function_, Named("method") = "L-BFGS-B",
        Named("data") = data_segment, Named("lower") = parameters_lower_bound_,
        Named("upper") = parameters_upper_bound_);
    result_coefficients_ = as<mat>(optim_result["par"]);
    result_residuals_ = mat();
    result_value_ = as<double>(optim_result["value"]);
  }
}

void Fastcpd::UpdateSenParametersStep(const int segment_start,
                                      const int segment_end, const int i) {
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

void Fastcpd::UpdateSenParametersSteps(const int segment_start,
                                       const unsigned int segment_end,
                                       const int i) {
  // This hack is to avoid the issue during momentum assigment of `solve`.
  colvec tmp = momentum_;
  const unsigned int multiple_epochs =
      as<int>((*multiple_epochs_function_)(segment_end - segment_start + 1));
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

double Fastcpd::UpdateCostValue(double value, const unsigned int nrows) {
  if (cost_adjustment_ == "MDL") {
    value = value * std::log2(M_E);
  }
  return value + GetCostAdjustmentValue(nrows);
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

void Fastcpd::UpdateSenParameters(const unsigned int t) {
  if (vanilla_percentage_ == 1) return;
  const int segment_index = index_max(find(segment_indices_ <= t - 1));
  colvec cum_coef_add = segment_coefficients_.row(segment_index).t(),
         coef_add = segment_coefficients_.row(segment_index).t();
  mat hessian_new;
  if (family_ == "binomial") {
    const rowvec x = data_.row(t - 1).tail(parameters_count_);
    const double prob = 1 / (1 + exp(-dot(coef_add, x)));
    hessian_new = (x.t() * x) * prob * (1 - prob);
  } else if (family_ == "poisson") {
    const rowvec x = data_.row(t - 1).tail(parameters_count_);
    cum_coef_add = coef_add;
    hessian_new = (x.t() * x) * exp(dot(coef_add, x));
  } else if (family_ == "lasso" || family_ == "gaussian") {
    const rowvec x = data_.row(t - 1).tail(parameters_count_);
    hessian_new = x.t() * x + epsilon_in_hessian_ * eye<mat>(parameters_count_,
                                                             parameters_count_);
  } else if (family_ == "arma") {
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

void Fastcpd::UpdateRClockTick(const std::string name) {
  if (!r_clock_.empty()) {
    rClock_.tick(name);
  }
}

void Fastcpd::UpdateRClockTock(const std::string name) {
  if (!r_clock_.empty()) {
    rClock_.tock(name);
  }
}

void Fastcpd::UpdateRProgress() {
  if (r_progress_) {
    rProgress_->tick();
  }
}

void Fastcpd::UpdateStep(unsigned int t) {
  UpdateRClockTick("pruning");
  UpdateSenParameters(t);
  colvec obj = GetObjectiveFunctionValues(t);

  // The following code is the manual implementation of `index_min` function
  // in Armadillo.
  objective_function_values_min_ = obj[0];
  objective_function_values_min_index_ = 0;
  for (unsigned int i = 1; i < pruned_set_size_; i++) {
    if (obj[i] < objective_function_values_min_) {
      objective_function_values_min_ = obj[i];
      objective_function_values_min_index_ = i;
    }
  }
  objective_function_values_(t) = objective_function_values_min_;
  change_points_[t] = pruned_set_[objective_function_values_min_index_];

  // The following code is the manual implementation of `find` function in
  // Armadillo.
  pruned_left_n_elem_ = 0;
  for (unsigned int i = 0; i < pruned_set_size_; i++) {
    if (obj[i] <=
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
  UpdateRClockTock("pruning");
  UpdateRProgress();
  checkUserInterrupt();
}

}  // namespace fastcpd::classes
