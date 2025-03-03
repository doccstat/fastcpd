#ifndef FASTCPD_CLASS_H_
#define FASTCPD_CLASS_H_

#include <memory>
#include <unordered_map>

// TODO(doccstat): Fix the order issue. It seems that `fastcpd_test.h` or
// `RcppArmadillo.h` should be included first.
#include "fastcpd_test.h"
#include "RProgress.h"
#include "RcppClock.h"

constexpr char kRProgress[] = "[:bar] :current/:total in :elapsed";

namespace fastcpd::classes {

class Fastcpd {
 public:
  Fastcpd(const double beta, const Rcpp::Nullable<Rcpp::Function> cost,
          const std::string cost_adjustment,
          const Rcpp::Nullable<Rcpp::Function> cost_gradient,
          const Rcpp::Nullable<Rcpp::Function> cost_hessian, const bool cp_only,
          const unsigned int d, const arma::mat data, const double epsilon,
          const std::string family,
          const Rcpp::Nullable<Rcpp::Function> multiple_epochs_function,
          const arma::colvec line_search, const arma::colvec lower,
          const double momentum_coef, const arma::colvec order, const int p,
          const unsigned int p_response, const double pruning_coef,
          const std::string r_clock, const bool r_progress,
          const int segment_count, const double trim, const arma::colvec upper,
          const double vanilla_percentage, const arma::mat variance_estimate,
          const bool warm_start);

  Rcpp::List Run();

 private:
  struct FunctionSet {
    arma::colvec (Fastcpd::*gradient)(const unsigned int segment_start,
                                      const unsigned int segment_end,
                                      const arma::colvec& theta);
    arma::mat (Fastcpd::*hessian)(const unsigned int segment_start,
                                  const unsigned int segment_end,
                                  const arma::colvec& theta);
    double (Fastcpd::*nll_sen)(const unsigned int segment_start,
                               const unsigned int segment_end,
                               arma::colvec theta);
    CostResult (Fastcpd::*nll_pelt)(const unsigned int segment_start,
                                    const unsigned int segment_end,
                                    const bool cv,
                                    const Rcpp::Nullable<arma::colvec>& start);
  };

  void CreateRClock(const std::string name);
  void CreateSenParameters();
  void CreateSegmentStatistics();
  void CreateSegmentStatisticsAndSenParameters();

  // Set \code{theta_sum} for a specific column.
  void CreateThetaSum(const unsigned int col, arma::colvec new_theta_sum);

  double GetCostAdjustmentValue(const unsigned nrows);

  // Solve logistic/poisson regression using Gradient Descent Extension to the
  // multivariate case
  // @return Negative log likelihood of the corresponding data with the given
  //   family.
  CostResult GetCostResult(const unsigned int segment_start,
                           const unsigned int segment_end,
                           Rcpp::Nullable<arma::colvec> theta,
                           const bool cv = false,
                           Rcpp::Nullable<arma::colvec> start = R_NilValue);
  Rcpp::List GetChangePointSet();
  double GetCostValue(const int tau, const unsigned int i, const int t);
  double GetCostValuePelt(const unsigned int segment_start,
                          const unsigned int segment_end, const unsigned int i);
  double GetCostValueSen(const unsigned int segment_start,
                         const unsigned int segment_end, const unsigned int i);
  arma::colvec GetGradientArma(const unsigned int segment_start,
                               const unsigned int segment_end,
                               const arma::colvec& theta);
  arma::colvec GetGradientBinomial(const unsigned int segment_start,
                                   const unsigned int segment_end,
                                   const arma::colvec& theta);
  arma::colvec GetGradientCustom(const unsigned int segment_start,
                                 const unsigned int segment_end,
                                 const arma::colvec& theta);
  arma::colvec GetGradientLm(const unsigned int segment_start,
                             const unsigned int segment_end,
                             const arma::colvec& theta);
  arma::colvec GetGradientMa(const unsigned int segment_start,
                             const unsigned int segment_end,
                             const arma::colvec& theta);
  arma::colvec GetGradientPoisson(const unsigned int segment_start,
                                  const unsigned int segment_end,
                                  const arma::colvec& theta);
  arma::mat GetHessianArma(const unsigned int segment_start,
                           const unsigned int segment_end,
                           const arma::colvec& theta);
  arma::mat GetHessianBinomial(const unsigned int segment_start,
                               const unsigned int segment_end,
                               const arma::colvec& theta);
  arma::mat GetHessianCustom(const unsigned int segment_start,
                             const unsigned int segment_end,
                             const arma::colvec& theta);
  arma::mat GetHessianLm(const unsigned int segment_start,
                         const unsigned int segment_end,
                         const arma::colvec& theta);
  arma::mat GetHessianMa(const unsigned int segment_start,
                         const unsigned int segment_end,
                         const arma::colvec& theta);
  arma::mat GetHessianPoisson(const unsigned int segment_start,
                              const unsigned int segment_end,
                              const arma::colvec& theta);
  CostResult GetNllPeltArma(const unsigned int segment_start,
                            const unsigned int segment_end, const bool cv,
                            const Rcpp::Nullable<arma::colvec>& start);
  CostResult GetNllPeltCustom(const unsigned int segment_start,
                              const unsigned int segment_end, const bool cv,
                              const Rcpp::Nullable<arma::colvec>& start);
  CostResult GetNllPeltGarch(const unsigned int segment_start,
                             const unsigned int segment_end, const bool cv,
                             const Rcpp::Nullable<arma::colvec>& start);
  CostResult GetNllPeltGlm(const unsigned int segment_start,
                           const unsigned int segment_end, const bool cv,
                           const Rcpp::Nullable<arma::colvec>& start);
  CostResult GetNllPeltLasso(const unsigned int segment_start,
                             const unsigned int segment_end, const bool cv,
                             const Rcpp::Nullable<arma::colvec>& start);
  CostResult GetNllPeltMean(const unsigned int segment_start,
                            const unsigned int segment_end, const bool cv,
                            const Rcpp::Nullable<arma::colvec>& start);
  CostResult GetNllPeltMeanVariance(const unsigned int segment_start,
                                    const unsigned int segment_end,
                                    const bool cv,
                                    const Rcpp::Nullable<arma::colvec>& start);
  CostResult GetNllPeltMgaussian(const unsigned int segment_start,
                                 const unsigned int segment_end, const bool cv,
                                 const Rcpp::Nullable<arma::colvec>& start);
  CostResult GetNllPeltVariance(const unsigned int segment_start,
                                const unsigned int segment_end, const bool cv,
                                const Rcpp::Nullable<arma::colvec>& start);
  double GetNllSenArma(const unsigned int segment_start,
                       const unsigned int segment_end, arma::colvec theta);
  double GetNllSenBinomial(const unsigned int segment_start,
                           const unsigned int segment_end, arma::colvec theta);
  double GetNllSenCustom(const unsigned int segment_start,
                         const unsigned int segment_end, arma::colvec theta);
  double GetNllSenLasso(const unsigned int segment_start,
                        const unsigned int segment_end, arma::colvec theta);
  double GetNllSenLm(const unsigned int segment_start,
                     const unsigned int segment_end, arma::colvec theta);
  double GetNllSenMa(const unsigned int segment_start,
                     const unsigned int segment_end, arma::colvec theta);
  double GetNllSenPoisson(const unsigned int segment_start,
                          const unsigned int segment_end, arma::colvec theta);
  arma::colvec GetObjectiveFunctionValues(unsigned int t);
  CostResult GetOptimizedCostResult(const unsigned int segment_start,
                                    const unsigned int segment_end);

  // Adjust cost value for MBIC and MDL.
  double UpdateCostValue(double value, const unsigned int nrows);
  arma::colvec UpdateChangePointSet();

  // Append new values to \code{fastcpd_parameters}.
  void UpdateSenParameters(const unsigned int t);
  void UpdateSenParameters(const unsigned int segment_start,
                           const unsigned int segment_end,
                           const unsigned int i);
  void UpdateSenParametersStep(const int segment_start, const int segment_end,
                               const int i);

  // Update the cost values for the segmentation.
  // @return A list containing new values of \code{theta_hat}, \code{theta_sum},
  //   \code{hessian}, and \code{momentum}.
  Rcpp::List UpdateSenParametersSteps(const int segment_start,
                                      const unsigned int segment_end,
                                      const int i, arma::colvec momentum);
  void UpdateStep(unsigned int t);

  // Append a new slice to \code{hessian}.
  void UpdateHessian(arma::mat new_hessian);

  // Prune the slices of \code{hessian}.
  void UpdateHessian(arma::ucolvec pruned_left);

  // Update \code{hessian} for a specific slice.
  void UpdateHessian(const unsigned int slice, arma::mat new_hessian);
  void UpdateMomentum(arma::colvec new_momentum);
  void UpdateRClockTick(const std::string name);
  void UpdateRClockTock(const std::string name);
  void UpdateRProgressStart();
  void UpdateRProgressTick();

  // Append a new column to \code{theta_hat}.
  void UpdateThetaHat(arma::colvec new_theta_hat);

  // Prune the columns of \code{theta_hat}.
  void UpdateThetaHat(arma::ucolvec pruned_left);

  // Update \code{theta_hat} for a specific column.
  void UpdateThetaHat(const unsigned int col, arma::colvec new_theta_hat);

  // Append a new column to \code{theta_sum}.
  void UpdateThetaSum(arma::colvec new_theta_sum);

  // Prune the columns of \code{theta_sum}.
  void UpdateThetaSum(arma::ucolvec pruned_left);

  // Update \code{theta_sum} for a specific column by adding to that column.
  void UpdateThetaSum(const unsigned int col, arma::colvec new_theta_sum);

  arma::colvec active_coefficients_count_;
  double beta_;
  arma::colvec change_points_;
  arma::mat coefficients_;
  arma::mat coefficients_sum_;
  const std::string cost_adjustment_;
  const std::unique_ptr<Rcpp::Function> cost_function_;
  const std::unique_ptr<Rcpp::Function> cost_gradient_;
  const std::unique_ptr<Rcpp::Function> cost_hessian_;
  const bool cp_only_;
  const arma::mat data_;
  const arma::mat data_c_;
  const unsigned int data_c_n_cols_;
  const unsigned int data_c_n_rows_;
  const double* data_c_ptr_;
  const unsigned int data_n_dims_;
  const unsigned int data_n_cols_;
  const unsigned int data_n_rows_;
  const double epsilon_in_hessian_;
  arma::colvec error_standard_deviation_;
  const std::string family_;
  static const std::unordered_map<std::string, FunctionSet>
      family_function_map_;
  arma::colvec (Fastcpd::*get_gradient_)(const unsigned int segment_start,
                                         const unsigned int segment_end,
                                         const arma::colvec& theta);
  arma::mat (Fastcpd::*get_hessian_)(const unsigned int segment_start,
                                     const unsigned int segment_end,
                                     const arma::colvec& theta);
  CostResult (Fastcpd::*get_nll_pelt_)(
      const unsigned int segment_start, const unsigned int segment_end,
      const bool cv, const Rcpp::Nullable<arma::colvec>& start);
  double (Fastcpd::*get_nll_sen_)(const unsigned int segment_start,
                                  const unsigned int segment_end,
                                  arma::colvec theta);
  arma::cube hessian_;
  double lasso_penalty_base_;
  arma::colvec line_search_;
  arma::colvec momentum_;
  const double momentum_coef_;
  const std::unique_ptr<Rcpp::Function> multiple_epochs_function_;
  arma::colvec objective_function_values_;
  double objective_function_values_min_;
  unsigned int objective_function_values_min_index_;
  const arma::colvec order_;
  const unsigned int parameters_count_;
  const arma::colvec parameters_lower_bound_;
  const arma::colvec parameters_upper_bound_;
  arma::ucolvec pruned_left_;
  unsigned int pruned_left_n_elem_;
  arma::ucolvec pruned_set_;
  unsigned int pruned_set_size_ = 2;
  const double pruning_coefficient_;
  const std::string r_clock_;
  const bool r_progress_;
  Rcpp::Clock rClock_;
  const unsigned int regression_response_count_;
  std::unique_ptr<RProgress::RProgress> rProgress_;
  arma::mat segment_coefficients_;
  const int segment_count_;
  arma::colvec segment_indices_;
  const double trim_;
  const bool use_warm_start_;
  const double vanilla_percentage_;
  const arma::mat variance_estimate_;
  arma::mat warm_start_;
  friend fastcpd::test::FastcpdTest;
};

}  // namespace fastcpd::classes

#endif  // FASTCPD_CLASS_H_
