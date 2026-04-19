#ifndef FASTCPD_TEMPLATE_H_
#define FASTCPD_TEMPLATE_H_

#ifndef NO_RCPP
#include <RcppArmadillo.h>
#include "RProgress.h"
#else
#include <armadillo>
#endif

#include "fastcpd_family.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace fastcpd {
namespace classes {

// Cost-function storage, conditionally compiled out for PELT-only families.
// SenFunctions<true>  — empty struct; empty-base optimisation gives it 0 bytes.
// SenFunctions<false> — holds all four std::function members (~160 bytes).
template <bool kIsPeltOnly>
struct SenFunctions {
  SenFunctions(
#ifndef NO_RCPP
      std::optional<Rcpp::Function> const& /* cost */,
#endif
      std::function<double(arma::mat)> const& /* cost_pelt */,
      std::function<double(arma::mat, arma::colvec)> const& /* cost_sen */,
      std::function<arma::colvec(arma::mat, arma::colvec)> const& /* cost_gradient */,
      std::function<arma::mat(arma::mat, arma::colvec)> const& /* cost_hessian */) {}
};

template <>
struct SenFunctions<false> {
#ifndef NO_RCPP
  SenFunctions(
      std::optional<Rcpp::Function> const& cost,
      std::function<double(arma::mat)> const& cost_pelt,
      std::function<double(arma::mat, arma::colvec)> const& cost_sen,
      std::function<arma::colvec(arma::mat, arma::colvec)> const& cost_gradient,
      std::function<arma::mat(arma::mat, arma::colvec)> const& cost_hessian)
      : cost_function_(cost),
        cost_function_pelt_(cost_pelt),
        cost_function_sen_(cost_sen),
        cost_gradient_(cost_gradient),
        cost_hessian_(cost_hessian) {}
  std::optional<Rcpp::Function> const cost_function_;
#else
  SenFunctions(
      std::function<double(arma::mat)> const& cost_pelt,
      std::function<double(arma::mat, arma::colvec)> const& cost_sen,
      std::function<arma::colvec(arma::mat, arma::colvec)> const& cost_gradient,
      std::function<arma::mat(arma::mat, arma::colvec)> const& cost_hessian)
      : cost_function_pelt_(cost_pelt),
        cost_function_sen_(cost_sen),
        cost_gradient_(cost_gradient),
        cost_hessian_(cost_hessian) {}
#endif
  std::function<double(arma::mat)> const cost_function_pelt_;
  std::function<double(arma::mat, arma::colvec)> const cost_function_sen_;
  std::function<arma::colvec(arma::mat, arma::colvec)> const cost_gradient_;
  std::function<arma::mat(arma::mat, arma::colvec)> const cost_hessian_;
};

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims = -1>
class Fastcpd : public SenFunctions<FamilyPolicy::is_pelt_only> {
 public:
  // Returns an empty matrix for PELT-only families (never used, avoids O(n)
  // allocation of coefficients/coefficients_sum for billion-row datasets).
  static arma::mat MakeSenMat(int const p, arma::uword const n) {
    if constexpr (FamilyPolicy::is_pelt_only) { return {}; }
    else { return arma::mat(p, n); }
  }
  // Same logic for the p×p×n Hessian cube.
  static arma::cube MakeHessianCube(int const p, arma::uword const n) {
    if constexpr (FamilyPolicy::is_pelt_only) { return {}; }
    else { return arma::cube(p, p, n); }
  }
  // warm_start_ is only touched inside if constexpr (supports_warm_start).
  static arma::mat MakeWarmStart(int const p, arma::uword const n) {
    if constexpr (!FamilyPolicy::supports_warm_start) { return {}; }
    else { return arma::zeros<arma::mat>(p, n); }
  }

  Fastcpd(
      double const beta,
#ifndef NO_RCPP
      std::optional<Rcpp::Function> const& cost,
#endif
      std::function<double(arma::mat)> const& cost_pelt,
      std::function<double(arma::mat, arma::colvec)> const& cost_sen,
      std::function<arma::colvec(arma::mat, arma::colvec)> const& cost_gradient,
      std::function<arma::mat(arma::mat, arma::colvec)> const& cost_hessian,
      bool const cp_only, arma::mat const& data, double const epsilon,
      std::string const& family,
      std::function<unsigned int(unsigned int)> const& multiple_epochs_function,
      arma::colvec const& line_search, arma::colvec const& lower,
      double const momentum_coef, arma::colvec const& order, int const p,
      unsigned int const p_response, double const pruning_coef,
      int const segment_count, double const trim,
      arma::colvec const& upper, double const vanilla_percentage,
      arma::mat const& variance_estimate, bool const warm_start)
      : SenFunctions<FamilyPolicy::is_pelt_only>(
#ifndef NO_RCPP
            cost,
#endif
            cost_pelt, cost_sen, cost_gradient, cost_hessian),
        active_coefficients_count_(arma::colvec(segment_count)),
        beta_(beta),
        change_points_(arma::zeros<arma::colvec>(data.n_rows + 1)),
        coefficients_(MakeSenMat(p, data.n_rows + 1)),
        coefficients_sum_(MakeSenMat(p, data.n_rows + 1)),
        cp_only_(cp_only),
        data_(data),
        data_c_(FamilyPolicy::CreateDataC(data, variance_estimate, p_response)),
        data_c_n_rows_(data_c_.n_rows),
        data_c_ptr_(data_c_.memptr()),
        data_n_dims_(FamilyPolicy::GetDataNDims(data)),
        data_n_rows_(data_.n_rows),
        epsilon_in_hessian_(epsilon),
        error_standard_deviation_(arma::colvec(segment_count)),
        family_(family),
        hessian_(MakeHessianCube(p, data.n_rows + 1)),
        line_search_(line_search),
        momentum_(arma::colvec(p)),
        momentum_coef_(momentum_coef),
        multiple_epochs_function_(multiple_epochs_function),
        objective_function_values_(arma::colvec(data_n_rows_ + 1)),
        objective_function_values_candidates_(arma::colvec(data_n_rows_ + 1)),
        objective_function_values_candidates_ptr_(
            objective_function_values_candidates_.memptr()),
        order_(order),
        parameters_count_(p),
        parameters_lower_bound_(lower),
        parameters_upper_bound_(upper),
        pruned_set_(arma::zeros<arma::ucolvec>(data_n_rows_ + 1)),
        pruning_coefficient_(pruning_coef),
        regression_response_count_(p_response),
#ifndef NO_RCPP
        rProgress_(kRProgress
            ? std::make_unique<RProgress::RProgress>(
                  "[:bar] :current/:total in :elapsed", data_n_rows_)
            : nullptr),
#endif
        segment_coefficients_(arma::mat(segment_count, p)),
        segment_count_(segment_count),
        segment_indices_(
            arma::round(arma::linspace(0, data_n_rows_, segment_count + 1))),
        trim_(trim),
        use_warm_start_(warm_start),
        vanilla_percentage_(vanilla_percentage),
        vanilla_end_(static_cast<unsigned int>(vanilla_percentage * data_n_rows_)),
        variance_estimate_(variance_estimate),
        variance_log_det_(arma::log_det_sympd(variance_estimate)),
        warm_start_(MakeWarmStart(p, data_n_rows_)) {
  }

  std::tuple<arma::colvec, arma::colvec, arma::colvec, arma::mat, arma::mat>
  Run();

  void CreateRProgress();
  void CreateSenParameters();
  void CreateSegmentStatistics();
  double GetCostAdjustmentValue(unsigned int const nrows);
  void GetCostResult(unsigned int const segment_start,
                     unsigned int const segment_end,
                     std::optional<arma::colvec> theta, bool const cv = false,
                     std::optional<arma::colvec> start = std::nullopt);
  std::tuple<arma::colvec, arma::colvec, arma::colvec, arma::mat, arma::mat>
  GetChangePointSet();
  double GetCostValue(int const tau, unsigned int const i);
  void GetCostValuePelt(unsigned int const segment_start,
                        unsigned int const segment_end, unsigned int const i);
  double GetCostValueSen(unsigned int const segment_start,
                         unsigned int const segment_end, unsigned int const i);

  void GetOptimizedCostResult(unsigned int const segment_start,
                              unsigned int const segment_end);
  arma::colvec UpdateChangePointSet();
  void UpdateSenParameters();
  void UpdateSenParametersStep(int const segment_start, int const segment_end,
                               int const i);
  void UpdateSenParametersSteps(int const segment_start,
                                unsigned int const segment_end, int const i);
  void UpdateStep();
  void UpdateRProgress();

  arma::colvec active_coefficients_count_;
  double beta_;
  arma::colvec change_points_;
  arma::mat coefficients_;
  arma::mat coefficients_sum_;
  bool const cp_only_;
  arma::mat const data_;
  arma::mat const data_c_;
  unsigned int const data_c_n_rows_;
  double const* data_c_ptr_;
  unsigned int const data_n_dims_;
  unsigned int const data_n_rows_;
  double const epsilon_in_hessian_;
  arma::colvec error_standard_deviation_;
  std::string const family_;
  arma::cube hessian_;
  double lasso_penalty_base_;
  arma::colvec line_search_;
  arma::colvec momentum_;
  double const momentum_coef_;
  std::function<unsigned int(unsigned int)> const multiple_epochs_function_;
  arma::colvec objective_function_values_;
  arma::colvec objective_function_values_candidates_;
  double* objective_function_values_candidates_ptr_;
  double objective_function_values_min_;
  unsigned int objective_function_values_min_index_;
  arma::colvec const order_;
  unsigned int const parameters_count_;
  arma::colvec const parameters_lower_bound_;
  arma::colvec const parameters_upper_bound_;
  unsigned int pruned_left_n_elem_;
  arma::ucolvec pruned_set_;
  unsigned int pruned_set_size_ = 2;
  double const pruning_coefficient_;
  unsigned int const regression_response_count_;
  arma::colvec result_coefficients_;
  arma::mat result_residuals_;
  double result_value_;
#ifndef NO_RCPP
  std::unique_ptr<RProgress::RProgress> rProgress_;
#endif
  arma::mat segment_coefficients_;
  int const segment_count_;
  arma::colvec segment_indices_;
  unsigned int t = 1;
  double const trim_;
  bool const use_warm_start_;
  double const vanilla_percentage_;
  // Precomputed integer boundary: avoids a float multiply per candidate in
  // GetCostValue when kVanillaOnly == false (SEN families only).
  unsigned int const vanilla_end_;
  arma::mat const variance_estimate_;
  // Precomputed once at construction — avoids calling log_det_sympd(Σ) once
  // per segment in mgaussian's GetNllPelt (hot path for VAR models).
  double const variance_log_det_;
  arma::mat warm_start_;
};

// --- Implementation ---

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
std::tuple<arma::colvec, arma::colvec, arma::colvec, arma::mat, arma::mat>
Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::Run() {
  CreateSegmentStatistics();
  CreateSenParameters();
  pruned_set_(1) = 1;
  objective_function_values_.fill(arma::datum::inf);
  objective_function_values_(0) = -beta_;
  objective_function_values_(1) = GetCostValue(0, 0);

  if constexpr (FamilyPolicy::has_optimized_run) {
    FamilyPolicy::Run(this);
  } else {
    CreateRProgress();
    UpdateRProgress();
    for (t = 2; t <= data_n_rows_; t++) {
      UpdateStep();
    }
  }
  return GetChangePointSet();
}

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
void Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::CreateRProgress() {
  if constexpr (kRProgress) {
#ifndef NO_RCPP
    rProgress_->tick(0);
#endif
  }
}

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
void Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::CreateSenParameters() {
  if constexpr (FamilyPolicy::is_pelt_only || kVanillaOnly) return;
  FamilyPolicy::CreateSenParameters(this);
}

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
void Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::CreateSegmentStatistics() {
  if constexpr (FamilyPolicy::has_optimized_run) return;
  if constexpr (FamilyPolicy::is_pelt_only) return;

  for (int segment_index = 0; segment_index < segment_count_; ++segment_index) {
    GetCostResult(segment_indices_(segment_index),
                  segment_indices_(segment_index + 1) - 1, std::nullopt, true,
                  std::nullopt);

    if (result_coefficients_.n_elem == (arma::uword)parameters_count_) {
      segment_coefficients_.row(segment_index) = result_coefficients_.t();
    }
    if constexpr (FamilyPolicy::has_error_std_dev) {
      arma::mat data_segment = data_.rows(segment_indices_(segment_index),
                                    segment_indices_(segment_index + 1) - 1);
      arma::colvec segment_residual =
          data_segment.col(0) -
          data_segment.cols(1, data_segment.n_cols - 1) * result_coefficients_;
      double err_var = arma::as_scalar(arma::mean(arma::square(segment_residual)));
      error_standard_deviation_(segment_index) = std::sqrt(err_var);
      active_coefficients_count_(segment_index) =
          arma::accu(arma::abs(result_coefficients_) > 0);
    }
  }
  lasso_penalty_base_ =
      arma::mean(error_standard_deviation_) * std::sqrt(2 * std::log(parameters_count_));
  if constexpr (FamilyPolicy::is_lasso) {
    beta_ = beta_ * (1 + arma::mean(active_coefficients_count_));
  }
}

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
double Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::GetCostAdjustmentValue(unsigned int const nrows) {
  if constexpr (kCostAdj == CostAdjustment::kBIC) {
    return 0;
  } else if constexpr (kCostAdj == CostAdjustment::kMBIC) {
    return parameters_count_ * std::log(nrows) / 2.0;
  } else {  // kMDL
    return parameters_count_ * std::log(nrows) / 2.0 * std::log2(M_E);
  }
}

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
void Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::GetCostResult(
    unsigned int const segment_start, unsigned int const segment_end,
    std::optional<arma::colvec> theta, bool const cv,
    std::optional<arma::colvec> start) {
  if (!theta.has_value()) {
    if constexpr (FamilyPolicy::is_pelt_only) {
      // PELT-only families: use fast path unless a full result is needed.
      if (cv) {
        FamilyPolicy::GetNllPelt(this, segment_start, segment_end, cv, start);
      } else {
        FamilyPolicy::GetNllPeltValue(this, segment_start, segment_end, cv, start);
      }
    } else {
      // kVanillaOnly=false guarantees vanilla_percentage_ < 1 at runtime,
      // so the full result is always needed; GetNllPeltValue is dead code here.
      if constexpr (kVanillaOnly) {
        if (cv) {
          FamilyPolicy::GetNllPelt(this, segment_start, segment_end, cv, start);
        } else {
          FamilyPolicy::GetNllPeltValue(this, segment_start, segment_end, cv, start);
        }
      } else {
        FamilyPolicy::GetNllPelt(this, segment_start, segment_end, cv, start);
      }
    }
  } else {
    result_coefficients_ = arma::colvec();
    result_residuals_ = arma::mat();
    arma::colvec const theta_ = theta.value();
    result_value_ = FamilyPolicy::GetNllSen(this, segment_start, segment_end, theta_);
  }
  if constexpr (kCostAdj == CostAdjustment::kMDL) {
    result_value_ = result_value_ * std::log2(M_E);
  }
  result_value_ += GetCostAdjustmentValue(segment_end - segment_start + 1);
}

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
std::tuple<arma::colvec, arma::colvec, arma::colvec, arma::mat, arma::mat>
Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::GetChangePointSet() {
  arma::colvec cp_set = UpdateChangePointSet();

  if (cp_only_) {
    return std::make_tuple(change_points_, cp_set, arma::colvec(), arma::mat(), arma::mat());
  }

  arma::colvec cp_loc_ = arma::zeros<arma::colvec>(cp_set.n_elem + 2);
  if (cp_set.n_elem) {
    cp_loc_.rows(1, cp_loc_.n_elem - 2) = cp_set;
  }
  cp_loc_(cp_loc_.n_elem - 1) = data_n_rows_;
  arma::colvec cp_loc = arma::unique(std::move(cp_loc_));
  arma::colvec cost_values = arma::zeros<arma::vec>(cp_loc.n_elem - 1);
  arma::mat thetas = arma::zeros<arma::mat>(parameters_count_, cp_loc.n_elem - 1);
  arma::mat residual = arma::zeros<arma::mat>(data_n_rows_, 1);
  unsigned int residual_next_start = 0;

  for (unsigned int i = 0; i < cp_loc.n_elem - 1; i++) {
    GetCostResult(cp_loc(i), cp_loc(i + 1) - 1, std::nullopt, true,
                  std::nullopt);
    cost_values(i) = result_value_;
    if constexpr (!FamilyPolicy::is_custom) {
      thetas.col(i) = result_coefficients_;
      if (result_residuals_.n_rows > 0) {
        if (result_residuals_.n_cols > 1 && residual.n_cols == 1) {
          residual = arma::zeros<arma::mat>(data_n_rows_, result_residuals_.n_cols);
        }
        residual.rows(residual_next_start,
                      residual_next_start + result_residuals_.n_rows - 1) =
            result_residuals_;
        residual_next_start += result_residuals_.n_rows;
      }
    }
  }
  return std::make_tuple(change_points_, cp_set, cost_values, residual, thetas);
}

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
double Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::GetCostValue(int const tau, unsigned int const i) {
  if constexpr (FamilyPolicy::is_pelt_only || kVanillaOnly) {
    GetCostValuePelt(tau, t - 1, i);
    return result_value_;
  } else {
    // vanilla_end_ is precomputed once at construction: integer compare vs
    // a float multiply per candidate per timestep.
    if (t > vanilla_end_) {
      return GetCostValueSen(tau, t - 1, i);
    } else {
      GetCostValuePelt(tau, t - 1, i);
      return result_value_;
    }
  }
}

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
void Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::GetCostValuePelt(
    unsigned int const segment_start, unsigned int const segment_end,
    unsigned int const i) {
  if constexpr (FamilyPolicy::supports_warm_start) {
    if (use_warm_start_ &&
        segment_end + 1 - segment_start >= 10 * parameters_count_) {
      GetCostResult(segment_start, segment_end, std::nullopt, false,
                    segment_coefficients_
                        .row(arma::index_max(arma::find(segment_indices_ <= segment_end)))
                        .t());
      warm_start_.col(segment_start) = result_coefficients_;
      if constexpr (!FamilyPolicy::is_pelt_only && !kVanillaOnly) {
        if (segment_end < vanilla_percentage_ * data_n_rows_) {
          coefficients_.col(i) = result_coefficients_;
          coefficients_sum_.col(i) += result_coefficients_;
        }
      }
      return;
    }
  }
  if constexpr (FamilyPolicy::is_pelt_only) {
    // Fast path for PELT-only families: bypass optional<colvec> construction
    // in GetCostResult so the compiler can inline the tight cumsum loop.
    // Receive the value as a return instead of a write through this-pointer;
    // the compiler keeps it in a register without a store-reload round trip.
    result_value_ = FamilyPolicy::template GetNllPeltValueFast<kNDims>(this, segment_start, segment_end);
    if constexpr (kCostAdj == CostAdjustment::kMDL) {
      result_value_ *= std::log2(M_E);
    }
    result_value_ += GetCostAdjustmentValue(segment_end - segment_start + 1);
  } else {
    GetCostResult(segment_start, segment_end, std::nullopt, false, std::nullopt);
    if constexpr (!kVanillaOnly) {
      if (segment_end < vanilla_percentage_ * data_n_rows_) {
        coefficients_.col(i) = result_coefficients_;
        coefficients_sum_.col(i) += result_coefficients_;
      }
    }
  }
}

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
double Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::GetCostValueSen(
    unsigned int const segment_start, unsigned int const segment_end,
    unsigned int const i) {
  unsigned int const segment_length = segment_end - segment_start + 1;
  double cval = 0;
  UpdateSenParametersSteps(segment_start, segment_end, i);
  arma::colvec theta = coefficients_sum_.col(i) / segment_length;
  if constexpr (FamilyPolicy::is_custom) {
    cval = FamilyPolicy::GetNllSen(this, segment_start, segment_end, theta);
  } else if constexpr (FamilyPolicy::is_lasso) {
    if (segment_length >= 3) {
      GetCostResult(segment_start, segment_end, theta, false, std::nullopt);
      cval = result_value_;
    }
  } else {
    if (segment_length >= parameters_count_) {
      GetCostResult(segment_start, segment_end, theta, false, std::nullopt);
      cval = result_value_;
    }
  }
  return cval;
}

#ifndef NO_RCPP
template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
void Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::GetOptimizedCostResult(
    unsigned int const segment_start, unsigned int const segment_end) {
  arma::mat const data_segment = data_.rows(segment_start, segment_end);
  if (parameters_count_ == 1) {
    Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
    Rcpp::Function optim = stats["optim"];
    Rcpp::List optim_result = optim(
        Rcpp::Named("par") = 0,
        Rcpp::Named("fn") = Rcpp::InternalFunction(
            +[](double theta, arma::mat data_, Rcpp::Function cost) {
              return cost(Rcpp::Named("data") = data_,
                          Rcpp::Named("theta") = std::log(theta / (1 - theta)));
            }),
        Rcpp::Named("method") = "Brent", Rcpp::Named("lower") = 0,
        Rcpp::Named("upper") = 1, Rcpp::Named("data") = data_segment,
        Rcpp::Named("cost") = this->cost_function_.value());
    arma::colvec par = Rcpp::as<arma::colvec>(optim_result["par"]);
    double value = Rcpp::as<double>(optim_result["value"]);
    result_coefficients_ = arma::log(par / (1 - par));
    result_residuals_ = arma::mat();
    result_value_ = std::exp(value) / (1 + std::exp(value));
  } else {
    Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
    Rcpp::Function optim = stats["optim"];
    Rcpp::List optim_result = optim(
        Rcpp::Named("par") = arma::zeros<arma::vec>(parameters_count_),
        Rcpp::Named("fn") = this->cost_function_.value(),
        Rcpp::Named("method") = "L-BFGS-B", Rcpp::Named("data") = data_segment,
        Rcpp::Named("lower") = parameters_lower_bound_,
        Rcpp::Named("upper") = parameters_upper_bound_);
    result_coefficients_ = Rcpp::as<arma::colvec>(optim_result["par"]);
    result_residuals_ = arma::mat();
    result_value_ = Rcpp::as<double>(optim_result["value"]);
  }
}
#else
template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
void Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::GetOptimizedCostResult(
    unsigned int const segment_start, unsigned int const segment_end) {
}
#endif

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
void Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::UpdateSenParametersStep(
    int const segment_start, int const segment_end, int const i) {
  arma::mat hessian_i = hessian_.slice(i);
  arma::colvec gradient;

  hessian_i +=
      FamilyPolicy::GetHessian(this, segment_start, segment_end, coefficients_.col(i));
  gradient =
      FamilyPolicy::GetGradient(this, segment_start, segment_end, coefficients_.col(i));

  arma::mat hessian_psd =
      hessian_i + epsilon_in_hessian_ *
                      arma::eye<arma::mat>(coefficients_.n_rows, coefficients_.n_rows);
  momentum_ = momentum_coef_ * momentum_ - arma::solve(hessian_psd, gradient);
  double best_learning_rate = 1;
  if constexpr (kLineSearch) {
    arma::colvec line_search_costs = arma::zeros<arma::colvec>(line_search_.n_elem);
    for (unsigned int line_search_index = 0;
         line_search_index < line_search_.n_elem; line_search_index++) {
      arma::colvec theta_candidate =
          coefficients_.col(i) + line_search_[line_search_index] * momentum_;
      arma::colvec theta_upper_bound =
          arma::min(std::move(theta_candidate), parameters_upper_bound_);
      arma::colvec theta_projected =
          arma::max(std::move(theta_upper_bound), parameters_lower_bound_);
      line_search_costs[line_search_index] =
          FamilyPolicy::GetNllSen(this, segment_start, segment_end, theta_projected);
    }
    best_learning_rate = line_search_[line_search_costs.index_min()];
  }
  coefficients_.col(i) += best_learning_rate * momentum_;
  coefficients_.col(i) =
      arma::min(coefficients_.col(i), parameters_upper_bound_);
  coefficients_.col(i) =
      arma::max(coefficients_.col(i), parameters_lower_bound_);

  if constexpr (FamilyPolicy::has_error_std_dev) {
    double hessian_norm = arma::norm(hessian_i, "fro");
    arma::vec normd = arma::abs(coefficients_.col(i));
    if constexpr (FamilyPolicy::is_lasso) {
      normd -= lasso_penalty_base_ / std::sqrt(segment_end - segment_start + 1) /
               hessian_norm;
    }
    coefficients_.col(i) = arma::sign(coefficients_.col(i)) %
                           arma::max(normd, arma::zeros<arma::colvec>(normd.n_elem));
  }

  hessian_.slice(i) = std::move(hessian_i);
}

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
void Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::UpdateSenParametersSteps(
    int const segment_start, unsigned int const segment_end, int const i) {
  arma::colvec tmp = momentum_;
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

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
arma::colvec Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::UpdateChangePointSet() {
  arma::colvec cp_set = arma::zeros<arma::colvec>(data_n_rows_);
  int ncpts = 0;
  int last = data_n_rows_;
  while (last != 0) {
    cp_set[ncpts] = last;
    last = change_points_[last];
    ncpts += 1;
  }
  cp_set = arma::sort(cp_set.rows(arma::find(cp_set > 0)));
  cp_set = cp_set(arma::find(cp_set > trim_ * data_n_rows_));
  cp_set = cp_set(arma::find(cp_set < (1 - trim_) * data_n_rows_));
  arma::colvec cp_set_ = arma::zeros<arma::vec>(cp_set.n_elem + 1);
  if (cp_set.n_elem) {
    cp_set_.rows(1, cp_set_.n_elem - 1) = std::move(cp_set);
  }
  cp_set = arma::sort(arma::unique(std::move(cp_set_)));

  arma::ucolvec cp_set_too_close = arma::find(arma::diff(cp_set) <= trim_ * data_n_rows_);
  if (cp_set_too_close.n_elem > 0) {
    int rest_element_count = cp_set.n_elem - cp_set_too_close.n_elem;
    arma::colvec cp_set_rest_left = arma::zeros<arma::vec>(rest_element_count),
           cp_set_rest_right = arma::zeros<arma::vec>(rest_element_count);
    for (unsigned int i = 0, i_left = 0, i_right = 0; i < cp_set.n_elem; i++) {
      if (arma::ucolvec left_find = arma::find(cp_set_too_close == i);
          left_find.n_elem == 0) {
        cp_set_rest_left(i_left) = cp_set(i);
        i_left++;
      }
      if (arma::ucolvec right_find = arma::find(cp_set_too_close == i - 1);
          right_find.n_elem == 0) {
        cp_set_rest_right(i_right) = cp_set(i);
        i_right++;
      }
    }
    cp_set = arma::floor((cp_set_rest_left + cp_set_rest_right) / 2);
  }
  return cp_set(arma::find(cp_set > 0));
}

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
void Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::UpdateSenParameters() {
  if constexpr (FamilyPolicy::is_pelt_only || kVanillaOnly) return;
  FamilyPolicy::UpdateSenParameters(this);
}

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
void Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::UpdateRProgress() {
  if constexpr (kRProgress) {
#ifndef NO_RCPP
    rProgress_->tick();
#endif
  }
}

template <typename FamilyPolicy, bool kRProgress, bool kVanillaOnly,
          CostAdjustment kCostAdj, bool kLineSearch, int kNDims>
void Fastcpd<FamilyPolicy, kRProgress, kVanillaOnly, kCostAdj, kLineSearch, kNDims>::UpdateStep() {
  UpdateSenParameters();

  // Cache raw pointers once — eliminates Armadillo bounds-checked operator()
  // on every candidate access in all four hot loops below.
  arma::uword const* const pruned_set_raw = pruned_set_.memptr();
  double const* const objfn_ptr = objective_function_values_.memptr();

  if constexpr (FamilyPolicy::is_pelt_only) {
    // PELT-only hot loop: split into a prefetch-issuing body and a short tail
    // so there is no per-iteration branch on the prefetch lookahead bound.
    // Memory latency ~200–400 cycles; ~20 cycles compute/candidate → distance
    // of 16 keeps the memory pipeline saturated.  Both the family's data_c_
    // slice and the objective function value are prefetched together.
    //
    // GetNllPeltValueFast + cost adjustment are inlined directly here to avoid
    // the GetCostValue → GetCostValuePelt → result_value_ (this-member) write
    // → GetCostValue read-back chain.  The compiler must treat writes through
    // `this` conservatively (potential aliasing), so the store-reload from
    // result_value_ cannot be eliminated even with full inlining.  Bypassing
    // it keeps the value in a register from GetNllPeltValueFast to the final
    // assignment with no round-trip through memory.
    // segment_length = (t-1) - s + 1 = t - s  (segment is [s, t-1]).
    auto eval_pelt_candidate = [this, objfn_ptr](unsigned int const s,
                                                  unsigned int const i) {
      double nll = FamilyPolicy::template GetNllPeltValueFast<kNDims>(this, s, t - 1);
      if constexpr (kCostAdj == CostAdjustment::kMBIC) {
        nll += static_cast<double>(parameters_count_) *
               std::log(static_cast<double>(t - s)) / 2.0;
      } else if constexpr (kCostAdj == CostAdjustment::kMDL) {
        nll = (nll + static_cast<double>(parameters_count_) *
                         std::log(static_cast<double>(t - s)) / 2.0) *
              std::log2(M_E);
      }
      objective_function_values_candidates_ptr_[i] = objfn_ptr[s] + nll + beta_;
    };
    constexpr unsigned int kPrefetchDist = 16;
    unsigned int const hot_end =
        pruned_set_size_ > kPrefetchDist ? pruned_set_size_ - kPrefetchDist : 0;
    // hot body — always prefetch, branch-free inner loop
    for (unsigned int i = 0; i < hot_end; i++) {
      unsigned int const next_s = pruned_set_raw[i + kPrefetchDist];
      FamilyPolicy::PrefetchCandidate(this, next_s);
#if defined(__GNUC__) || defined(__clang__)
      __builtin_prefetch(objfn_ptr + next_s, 0, 0);
#endif
      eval_pelt_candidate(pruned_set_raw[i], i);
    }
    // tail — last kPrefetchDist candidates, data already en route or in cache
    for (unsigned int i = hot_end; i < pruned_set_size_; i++) {
      eval_pelt_candidate(pruned_set_raw[i], i);
    }
  } else {
    // SEN family candidate loop.
    // When !kVanillaOnly the last candidate always skips GetCostValue (it
    // represents the trivially-bounded new segment).  Peel that final
    // iteration out of the loop body so there is no per-iteration branch.
    // Note: i + 1 < size avoids unsigned underflow when size == 0.
    if constexpr (kVanillaOnly) {
      for (unsigned int i = 0; i < pruned_set_size_; i++) {
        unsigned int const s = pruned_set_raw[i];
        objective_function_values_candidates_ptr_[i] =
            objfn_ptr[s] + GetCostValue(s, i) + beta_;
      }
    } else {
      for (unsigned int i = 0; i + 1 < pruned_set_size_; i++) {
        unsigned int const s = pruned_set_raw[i];
        objective_function_values_candidates_ptr_[i] =
            objfn_ptr[s] + GetCostValue(s, i) + beta_;
      }
      objective_function_values_candidates_ptr_[pruned_set_size_ - 1] =
          objfn_ptr[pruned_set_raw[pruned_set_size_ - 1]] + beta_;
    }
  }

  // SIMD-friendly argmin: std::min_element auto-vectorizes (no data-dependent
  // branch per iteration), then a single pointer subtraction gives the index.
  double const* const min_it = std::min_element(
      objective_function_values_candidates_ptr_,
      objective_function_values_candidates_ptr_ + pruned_set_size_);
  objective_function_values_min_ = *min_it;
  objective_function_values_min_index_ =
      static_cast<unsigned int>(min_it - objective_function_values_candidates_ptr_);
  objective_function_values_(t) = objective_function_values_min_;
  change_points_[t] = pruned_set_raw[objective_function_values_min_index_];

  // Precompute threshold once — eliminates a float add+subtract from every
  // iteration of the pruning loop (was recomputed per candidate before).
  double const pruning_threshold =
      objective_function_values_min_ + beta_ - pruning_coefficient_;
  pruned_left_n_elem_ = 0;
  for (unsigned int i = 0; i < pruned_set_size_; i++) {
    if (objective_function_values_candidates_ptr_[i] <= pruning_threshold) {
      pruned_set_[pruned_left_n_elem_] = pruned_set_raw[i];
      if constexpr (!FamilyPolicy::is_pelt_only && !kVanillaOnly) {
        if (pruned_left_n_elem_ != i) {
          std::memcpy(coefficients_.colptr(pruned_left_n_elem_),
                 coefficients_.colptr(i), sizeof(double) * parameters_count_);
          std::memcpy(coefficients_sum_.colptr(pruned_left_n_elem_),
                 coefficients_sum_.colptr(i), sizeof(double) * parameters_count_);
          std::memcpy(hessian_.slice(pruned_left_n_elem_).memptr(),
                 hessian_.slice(i).memptr(),
                 sizeof(double) * parameters_count_ * parameters_count_);
        }
      }
      pruned_left_n_elem_++;
    }
  }
  pruned_set_size_ = pruned_left_n_elem_;
  pruned_set_[pruned_set_size_] = t;
  pruned_set_size_++;
  UpdateRProgress();
}

} // namespace classes
} // namespace fastcpd

#endif
