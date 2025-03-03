#include "fastcpd_class.h"

#include "RProgress.h"

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
using ::std::string;
using ::std::string_view;
using ::std::unique_ptr;
using ::std::unordered_map;
using ::std::vector;

using ::arma::vec;
using ::Rcpp::Environment;
using ::Rcpp::InternalFunction;

using ::Rcpp::Rcout;

#define ERROR(msg) \
  Rcout << "error: " << __FILE__ << ": " << __LINE__ << ": " << msg << std::endl
#define FATAL(msg)                                                  \
  Rcout << "fatal: " << __FILE__ << ": " << __LINE__ << ": " << msg \
        << std::endl;                                               \
  throw std::runtime_error(msg)

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
      coefficients_(mat(p, 1)),
      coefficients_sum_(mat(p, 1)),
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
      hessian_(cube(p, p, 1)),
      line_search_(line_search),
      momentum_(vec(p)),
      momentum_coef_(momentum_coef),
      multiple_epochs_function_([&]() -> unique_ptr<Function> {
        if (multiple_epochs_function.isNotNull()) {
          return make_unique<Function>(multiple_epochs_function);
        }
        return nullptr;
      }()),
      order_(order),
      parameters_count_(p),
      parameters_lower_bound_(lower),
      parameters_upper_bound_(upper),
      pruned_left_(ucolvec(data_n_rows_ + 1)),
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
      const FunctionSet& func_set = it->second;
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
  ucolvec r_t_set = zeros<ucolvec>(data_n_rows_ + 1);
  r_t_set(1) = 1;
  unsigned int r_t_count = 2;

  colvec cp_sets = zeros<colvec>(data_n_rows_ + 1);
  colvec fvec = zeros<vec>(data_n_rows_ + 1);
  fvec.fill(arma::datum::inf);
  fvec(0) = -beta_;
  fvec(1) = GetCostValue(0, 0, 1);

  if (family_ == "mean" && cost_adjustment_ == "MBIC") {
    double* obj = (double*)calloc(data_n_rows_ + 1, sizeof(double));
    double two_norm;
    unsigned int i, pi;

    for (unsigned int t = 2; t <= data_n_rows_; t++) {
      for (i = 0; i < r_t_count; i++) {
        two_norm = (data_c_ptr_[t] - data_c_ptr_[r_t_set[i]]) *
                   (data_c_ptr_[t] - data_c_ptr_[r_t_set[i]]);
        for (pi = 1; pi < parameters_count_; pi++) {
          two_norm += (data_c_ptr_[t + data_c_n_rows_ * pi] -
                       data_c_ptr_[r_t_set[i] + data_c_n_rows_ * pi]) *
                      (data_c_ptr_[t + data_c_n_rows_ * pi] -
                       data_c_ptr_[r_t_set[i] + data_c_n_rows_ * pi]);
        }
        obj[i] =
            fvec[r_t_set[i]] +
            ((data_c_ptr_[t + data_c_n_rows_ * parameters_count_] -
              data_c_ptr_[r_t_set[i] + data_c_n_rows_ * parameters_count_]) -
             two_norm / (t - r_t_set[i])) /
                2.0 +
            std::log(t - r_t_set[i]) / 2.0 + beta_;
      }

      min_objective_function_value_ = obj[0];
      min_objective_function_value_index_ = 0;
      for (i = 1; i < r_t_count; i++) {
        if (obj[i] < min_objective_function_value_) {
          min_objective_function_value_ = obj[i];
          min_objective_function_value_index_ = i;
        }
      }
      fvec(t) = min_objective_function_value_;
      cp_sets[t] = r_t_set[min_objective_function_value_index_];

      pruned_left_n_elem_ = 0;
      for (i = 0; i < r_t_count; i++) {
        if (obj[i] <=
            min_objective_function_value_ + beta_ - pruning_coefficient_) {
          r_t_set[pruned_left_n_elem_] = r_t_set[i];
          pruned_left_n_elem_++;
        }
      }
      r_t_count = pruned_left_n_elem_;
      r_t_set[r_t_count] = t;
      r_t_count++;
    }
  } else {
    CreateSegmentStatisticsAndSenParameters();
    for (unsigned int t = 2; t <= data_n_rows_; t++) {
      UpdateStep(t, r_t_set, r_t_count, cp_sets, fvec);
    }
  }

  List result = GetChangePointSet(cp_sets);

  return result;
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
    rowvec segment_theta =
        GetCostResult(segment_indices_(segment_index),
                      segment_indices_(segment_index + 1) - 1, R_NilValue, true,
                      R_NilValue)
            .par;

    // Initialize the estimated coefficients for each segment to be the
    // estimated coefficients in the segment.
    segment_coefficients_.row(segment_index) = segment_theta;
    if (family_ == "lasso" || family_ == "gaussian") {
      mat data_segment = data_.rows(segment_indices_(segment_index),
                                    segment_indices_(segment_index + 1) - 1);
      colvec segment_residual =
          data_segment.col(0) -
          data_segment.cols(1, data_segment.n_cols - 1) * segment_theta.t();
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

void Fastcpd::CreateSegmentStatisticsAndSenParameters() {
  CreateSegmentStatistics();
  CreateSenParameters();
  checkUserInterrupt();
  UpdateRProgressStart();
  UpdateRProgressTick();
}

void Fastcpd::CreateThetaSum(const unsigned int col, colvec new_theta_sum) {
  coefficients_sum_.col(col) = new_theta_sum;
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

CostResult Fastcpd::GetCostResult(const unsigned int segment_start,
                                  const unsigned int segment_end,
                                  Nullable<colvec> theta, const bool cv,
                                  Nullable<colvec> start) {
  CostResult cost_result;
  if (theta.isNull()) {
    cost_result = (this->*get_nll_pelt_)(segment_start, segment_end, cv, start);
  } else {
    cost_result = CostResult{
        {colvec()},
        {colvec()},
        (this->*get_nll_sen_)(segment_start, segment_end, as<colvec>(theta))};
  }
  cost_result.value =
      UpdateCostValue(cost_result.value, segment_end - segment_start + 1);
  return cost_result;
}

List Fastcpd::GetChangePointSet(const colvec raw_cp_set) {
  CreateRClock(r_clock_);

  colvec cp_set = UpdateChangePointSet(raw_cp_set);

  if (cp_only_) {
    return List::create(
        Named("raw_cp_set") = raw_cp_set, Named("cp_set") = cp_set,
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
    CostResult cost_result = GetCostResult(cp_loc(i), cp_loc(i + 1) - 1,
                                           R_NilValue, false, R_NilValue);

    cost_values(i) = cost_result.value;

    // Parameters are not involved for PELT.
    if (vanilla_percentage_ < 1 || family_ == "garch") {
      thetas.col(i) = colvec(cost_result.par);
    }

    // Residual is only calculated for built-in families.
    if (family_ != "custom" && family_ != "mean" && family_ != "variance" &&
        family_ != "meanvariance") {
      mat cost_optim_residual = cost_result.residuals;
      residual.rows(residual_next_start,
                    residual_next_start + cost_optim_residual.n_rows - 1) =
          cost_optim_residual;
      residual_next_start += cost_optim_residual.n_rows;
    }
  }
  return List::create(Named("raw_cp_set") = raw_cp_set,
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
  CostResult cost_result;
  if ((family_ == "binomial" || family_ == "poisson") &&
      (use_warm_start_ &&
       segment_end + 1 - segment_start >= 10 * parameters_count_)) {
    cost_result = GetCostResult(
        segment_start, segment_end, R_NilValue, false,
        wrap(segment_coefficients_
                 .row(index_max(find(segment_indices_ <= segment_end)))
                 .t())
        // Or use `wrap(start.col(segment_start))` for warm start.
    );
    warm_start_.col(segment_start) = colvec(cost_result.par);
  } else {
    cost_result = GetCostResult(segment_start, segment_end, R_NilValue, false,
                                R_NilValue);
  }
  cval = cost_result.value;

  // If `vanilla_percentage_` is not 1, then we need to keep track of
  // thetas for later `fastcpd` steps.
  if (vanilla_percentage_ < 1 &&
      segment_end < vanilla_percentage_ * data_n_rows_) {
    UpdateThetaHat(i, cost_result.par);
    UpdateThetaSum(i, cost_result.par);
  }
  return cval;
}

double Fastcpd::GetCostValueSen(const unsigned int segment_start,
                                const unsigned int segment_end,
                                const unsigned int i) {
  const unsigned int segment_length = segment_end - segment_start + 1;
  double cval = 0;
  UpdateSenParameters(segment_start, segment_end, i);
  colvec theta = coefficients_sum_.col(i) / segment_length;
  if (family_ == "custom") {
    cval = (this->*get_nll_sen_)(segment_start, segment_end, theta);
  } else if ((family_ != "lasso" && segment_length >= parameters_count_) ||
             (family_ == "lasso" && segment_length >= 3)) {
    cval = GetCostResult(segment_start, segment_end, wrap(theta), false,
                         R_NilValue)
               .value;
  }
  // else segment_length < parameters_count_ or for lasso segment_length < 3
  return cval;
}

colvec Fastcpd::GetObjectiveFunctionValues(const colvec& fvec,
                                           const ucolvec& r_t_set,
                                           unsigned int r_t_count,
                                           unsigned int t) {
  colvec cval = zeros<vec>(r_t_count);
  UpdateRClockTick("r_t_set_for_loop");
  unsigned int loop_end = r_t_count - (vanilla_percentage_ != 1);
  for (unsigned int i = 0; i < loop_end; i++) {
    cval(i) = GetCostValue(r_t_set(i), i, t);
  }
  UpdateRClockTock("r_t_set_for_loop");
  colvec obj = cval + fvec.rows(r_t_set.rows(0, r_t_count - 1)) + beta_;
  return obj;
}

CostResult Fastcpd::GetOptimizedCostResult(const unsigned int segment_start,
                                           const unsigned int segment_end) {
  CostResult cost_result;
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
    cost_result = {
        {log(par / (1 - par))}, {colvec()}, exp(value) / (1 + exp(value))};
  } else {
    Environment stats = Environment::namespace_env("stats");
    Function optim = stats["optim"];
    List optim_result = optim(
        Named("par") = zeros<vec>(parameters_count_),
        Named("fn") = *cost_function_, Named("method") = "L-BFGS-B",
        Named("data") = data_segment, Named("lower") = parameters_lower_bound_,
        Named("upper") = parameters_upper_bound_);
    cost_result = {
        {as<colvec>(optim_result["par"])}, {colvec()}, optim_result["value"]};
  }
  return cost_result;
}

void Fastcpd::UpdateSenParameters(const unsigned int segment_start,
                                  const unsigned int segment_end,
                                  const unsigned int i) {
  List cost_update_result =
      UpdateSenParametersSteps(segment_start, segment_end, i, momentum_);
  UpdateThetaHat(i, as<colvec>(cost_update_result[0]));
  CreateThetaSum(i, as<colvec>(cost_update_result[1]));
  UpdateHessian(i, as<mat>(cost_update_result[2]));
  UpdateMomentum(as<colvec>(cost_update_result[3]));
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

  // Calculate momentum step
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

  // Update coefficients_ with momentum
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

List Fastcpd::UpdateSenParametersSteps(const int segment_start,
                                       const unsigned int segment_end,
                                       const int i, colvec momentum) {
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
  return List::create(coefficients_.col(i), coefficients_sum_.col(i),
                      hessian_.slice(i), momentum);
}

double Fastcpd::UpdateCostValue(double value, const unsigned int nrows) {
  if (cost_adjustment_ == "MDL") {
    value = value * std::log2(M_E);
  }
  return value + GetCostAdjustmentValue(nrows);
}

colvec Fastcpd::UpdateChangePointSet(const colvec raw_cp_set) {
  // Remove change points close to the boundaries.
  colvec cp_set = zeros<colvec>(data_n_rows_);
  int ncpts = 0;
  int last = data_n_rows_;
  while (last != 0) {
    cp_set[ncpts] = last;
    last = raw_cp_set[last];
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
  const int segment_index = index_max(find(segment_indices_ <= t - 1));
  rowvec cum_coef_add = segment_coefficients_.row(segment_index),
         coef_add = segment_coefficients_.row(segment_index);
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
    hessian_new = (this->*get_hessian_)(0, t - 1, coef_add.t());
  } else if (family_ == "custom") {
    hessian_new = zeros<mat>(parameters_count_, parameters_count_);
  }

  UpdateThetaHat(coef_add.t());
  UpdateThetaSum(cum_coef_add.t());
  UpdateHessian(hessian_new);
}

void Fastcpd::UpdateHessian(const unsigned int slice, mat new_hessian) {
  hessian_.slice(slice) = new_hessian;
}

void Fastcpd::UpdateHessian(mat new_hessian) {
  hessian_ = join_slices(hessian_, new_hessian);
}

void Fastcpd::UpdateHessian(ucolvec pruned_left) {
  hessian_ = hessian_.slices(pruned_left);
}

void Fastcpd::UpdateMomentum(colvec new_momentum) { momentum_ = new_momentum; }

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

void Fastcpd::UpdateRProgressStart() {
  if (r_progress_) {
    rProgress_->tick(0);
  }
}

void Fastcpd::UpdateRProgressTick() {
  if (r_progress_) {
    rProgress_->tick();
  }
}

void Fastcpd::UpdateStep(unsigned int t, ucolvec& r_t_set,
                         unsigned int& r_t_count, colvec& cp_sets,
                         colvec& fvec) {
  UpdateRClockTick("pruning");

  if (vanilla_percentage_ != 1) {
    UpdateSenParameters(t);
  }

  colvec obj = GetObjectiveFunctionValues(fvec, r_t_set, r_t_count, t);

  // The following code is the manual implementation of `index_min` function
  // in Armadillo.
  min_objective_function_value_ = obj[0];
  min_objective_function_value_index_ = 0;
  for (unsigned int i = 1; i < r_t_count; i++) {
    if (obj[i] < min_objective_function_value_) {
      min_objective_function_value_ = obj[i];
      min_objective_function_value_index_ = i;
    }
  }
  fvec(t) = min_objective_function_value_;
  cp_sets[t] = r_t_set[min_objective_function_value_index_];

  // The following code is the manual implementation of `find` function in
  // Armadillo.
  pruned_left_n_elem_ = 0;
  for (unsigned int i = 0; i < r_t_count; i++) {
    if (obj[i] <=
        min_objective_function_value_ + beta_ - pruning_coefficient_) {
      r_t_set[pruned_left_n_elem_] = r_t_set[i];
      pruned_left_[pruned_left_n_elem_] = i;
      pruned_left_n_elem_++;
    }
  }
  r_t_count = pruned_left_n_elem_;
  r_t_set[r_t_count] = t;
  r_t_count++;

  if (vanilla_percentage_ != 1) {
    UpdateThetaHat(pruned_left_.rows(0, pruned_left_n_elem_ - 1));
    UpdateThetaSum(pruned_left_.rows(0, pruned_left_n_elem_ - 1));
    UpdateHessian(pruned_left_.rows(0, pruned_left_n_elem_ - 1));
  }

  UpdateRClockTock("pruning");

  checkUserInterrupt();
  UpdateRProgressTick();
}

void Fastcpd::UpdateThetaHat(colvec new_theta_hat) {
  coefficients_ = join_rows(coefficients_, new_theta_hat);
}

void Fastcpd::UpdateThetaHat(const unsigned int col, colvec new_theta_hat) {
  coefficients_.col(col) = new_theta_hat;
}

void Fastcpd::UpdateThetaHat(ucolvec pruned_left) {
  coefficients_ = coefficients_.cols(pruned_left);
}

void Fastcpd::UpdateThetaSum(colvec new_theta_sum) {
  coefficients_sum_ = join_rows(coefficients_sum_, new_theta_sum);
}

void Fastcpd::UpdateThetaSum(ucolvec pruned_left) {
  coefficients_sum_ = coefficients_sum_.cols(pruned_left);
}

void Fastcpd::UpdateThetaSum(const unsigned int col, colvec new_theta_sum) {
  coefficients_sum_.col(col) += new_theta_sum;
}

}  // namespace fastcpd::classes
