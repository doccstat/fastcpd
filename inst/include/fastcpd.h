#ifndef FASTCPD_CLASS_H_
#define FASTCPD_CLASS_H_

#include <RcppArmadillo.h>

#include <memory>
#include <unordered_map>

#include "RProgress.h"

#define ERROR(msg)                                                        \
  Rcpp::Rcout << "error: " << __FILE__ << ": " << __LINE__ << ": " << msg \
              << std::endl
#define FATAL(msg)                                                        \
  Rcpp::Rcout << "fatal: " << __FILE__ << ": " << __LINE__ << ": " << msg \
              << std::endl;                                               \
  throw std::runtime_error(msg)

constexpr char kRProgress[] = "[:bar] :current/:total in :elapsed";

namespace fastcpd {

namespace test {

// The FastcpdTest class provides static helper methods for testing Fastcpd.
class FastcpdTest {
 public:
  // Computes the gradient using the Armadillo implementation.
  static arma::colvec GetGradientArma(const arma::mat& data,
                                      const unsigned int segment_start,
                                      const unsigned int segment_end,
                                      const arma::colvec& theta);

  // Computes the Hessian using the Armadillo implementation.
  static arma::mat GetHessianArma(const arma::mat& data,
                                  const unsigned int segment_start,
                                  const unsigned int segment_end,
                                  const arma::colvec& theta);

  // Computes the Hessian for the binomial model.
  static arma::mat GetHessianBinomial(const arma::mat& data,
                                      const unsigned int segment_start,
                                      const unsigned int segment_end,
                                      const arma::colvec& theta);

  // Computes the Hessian for the Poisson model.
  static arma::mat GetHessianPoisson(const arma::mat& data,
                                     const unsigned int segment_start,
                                     const unsigned int segment_end,
                                     const arma::colvec& theta);

  // Computes the negative log-likelihood for the SEN model.
  static double GetNllSen(const arma::mat& data,
                          const unsigned int segment_start,
                          const unsigned int segment_end,
                          const arma::colvec& theta);

  // Computes the negative log-likelihood for the PELT model.
  static std::tuple<arma::mat, arma::mat, double> GetNllPelt(
      const arma::mat& data, const unsigned int segment_start,
      const unsigned int segment_end, const bool cv,
      const Rcpp::Nullable<arma::colvec>& start);
};

}  // namespace test

namespace classes {

class Fastcpd {
 public:
  Fastcpd(const double beta, const Rcpp::Nullable<Rcpp::Function>& cost,
          const std::string& cost_adjustment,
          const Rcpp::Nullable<Rcpp::Function>& cost_gradient,
          const Rcpp::Nullable<Rcpp::Function>& cost_hessian,
          const bool cp_only, const arma::mat& data, const double epsilon,
          const std::string& family,
          const Rcpp::Nullable<Rcpp::Function>& multiple_epochs_function,
          const arma::colvec& line_search, const arma::colvec& lower,
          const double momentum_coef, const arma::colvec& order, const int p,
          const unsigned int p_response, const double pruning_coef,
          const bool r_progress, const int segment_count, const double trim,
          const arma::colvec& upper, const double vanilla_percentage,
          const arma::mat& variance_estimate, const bool warm_start);

  Rcpp::List Run();

 private:
  struct FunctionSet {
    arma::colvec (Fastcpd::*gradient)(const unsigned int segment_start,
                                      const unsigned int segment_end,
                                      const arma::colvec& theta);
    arma::mat (Fastcpd::*hessian)(const unsigned int segment_start,
                                  const unsigned int segment_end,
                                  const arma::colvec& theta);
    void (Fastcpd::*nll_pelt)(const unsigned int segment_start,
                              const unsigned int segment_end, const bool cv,
                              const Rcpp::Nullable<arma::colvec>& start);
    void (Fastcpd::*nll_pelt_value)(const unsigned int segment_start,
                                    const unsigned int segment_end,
                                    const bool cv,
                                    const Rcpp::Nullable<arma::colvec>& start);
    double (Fastcpd::*nll_sen)(const unsigned int segment_start,
                               const unsigned int segment_end,
                               const arma::colvec& theta);
  };

  void CreateRProgress();
  void CreateSenParameters();
  void CreateSegmentStatistics();
  double GetCostAdjustmentValue(const unsigned int nrows);
  void GetCostResult(const unsigned int segment_start,
                     const unsigned int segment_end,
                     Rcpp::Nullable<arma::colvec> theta, const bool cv = false,
                     Rcpp::Nullable<arma::colvec> start = R_NilValue);
  Rcpp::List GetChangePointSet();
  double GetCostValue(const int tau, const unsigned int i);
  void GetCostValuePelt(const unsigned int segment_start,
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
  void GetNllPeltArma(const unsigned int segment_start,
                      const unsigned int segment_end, const bool cv,
                      const Rcpp::Nullable<arma::colvec>& start);
  void GetNllPeltCustom(const unsigned int segment_start,
                        const unsigned int segment_end, const bool cv,
                        const Rcpp::Nullable<arma::colvec>& start);
  void GetNllPeltGarch(const unsigned int segment_start,
                       const unsigned int segment_end, const bool cv,
                       const Rcpp::Nullable<arma::colvec>& start);
  void GetNllPeltGlm(const unsigned int segment_start,
                     const unsigned int segment_end, const bool cv,
                     const Rcpp::Nullable<arma::colvec>& start);
  void GetNllPeltLasso(const unsigned int segment_start,
                       const unsigned int segment_end, const bool cv,
                       const Rcpp::Nullable<arma::colvec>& start);
  void GetNllPeltMean(const unsigned int segment_start,
                      const unsigned int segment_end, const bool cv,
                      const Rcpp::Nullable<arma::colvec>& start);
  void GetNllPeltMeanValue(const unsigned int segment_start,
                           const unsigned int segment_end, const bool cv,
                           const Rcpp::Nullable<arma::colvec>& start);
  void GetNllPeltMeanvariance(const unsigned int segment_start,
                              const unsigned int segment_end, const bool cv,
                              const Rcpp::Nullable<arma::colvec>& start);
  void GetNllPeltMeanvarianceValue(const unsigned int segment_start,
                                   const unsigned int segment_end,
                                   const bool cv,
                                   const Rcpp::Nullable<arma::colvec>& start);
  void GetNllPeltMgaussian(const unsigned int segment_start,
                           const unsigned int segment_end, const bool cv,
                           const Rcpp::Nullable<arma::colvec>& start);
  void GetNllPeltVariance(const unsigned int segment_start,
                          const unsigned int segment_end, const bool cv,
                          const Rcpp::Nullable<arma::colvec>& start);
  void GetNllPeltVarianceValue(const unsigned int segment_start,
                               const unsigned int segment_end, const bool cv,
                               const Rcpp::Nullable<arma::colvec>& start);
  double GetNllSenArma(const unsigned int segment_start,
                       const unsigned int segment_end,
                       const arma::colvec& theta);
  double GetNllSenBinomial(const unsigned int segment_start,
                           const unsigned int segment_end,
                           const arma::colvec& theta);
  double GetNllSenCustom(const unsigned int segment_start,
                         const unsigned int segment_end,
                         const arma::colvec& theta);
  double GetNllSenLasso(const unsigned int segment_start,
                        const unsigned int segment_end,
                        const arma::colvec& theta);
  double GetNllSenLm(const unsigned int segment_start,
                     const unsigned int segment_end, const arma::colvec& theta);
  double GetNllSenMa(const unsigned int segment_start,
                     const unsigned int segment_end, const arma::colvec& theta);
  double GetNllSenPoisson(const unsigned int segment_start,
                          const unsigned int segment_end,
                          const arma::colvec& theta);
  void GetOptimizedCostResult(const unsigned int segment_start,
                              const unsigned int segment_end);
  arma::colvec UpdateChangePointSet();
  void UpdateSenParameters();
  void UpdateSenParametersStep(const int segment_start, const int segment_end,
                               const int i);
  void UpdateSenParametersSteps(const int segment_start,
                                const unsigned int segment_end, const int i);
  void UpdateStep();
  void UpdateRProgress();

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
  arma::colvec (Fastcpd::* const get_gradient_)(
      const unsigned int segment_start, const unsigned int segment_end,
      const arma::colvec& theta);
  arma::mat (Fastcpd::* const get_hessian_)(const unsigned int segment_start,
                                            const unsigned int segment_end,
                                            const arma::colvec& theta);
  void (Fastcpd::* const get_nll_pelt_)(
      const unsigned int segment_start, const unsigned int segment_end,
      const bool cv, const Rcpp::Nullable<arma::colvec>& start);
  void (Fastcpd::* const get_nll_pelt_value_)(
      const unsigned int segment_start, const unsigned int segment_end,
      const bool cv, const Rcpp::Nullable<arma::colvec>& start);
  double (Fastcpd::* const get_nll_sen_)(const unsigned int segment_start,
                                         const unsigned int segment_end,
                                         const arma::colvec& theta);
  arma::cube hessian_;
  double lasso_penalty_base_;
  arma::colvec line_search_;
  arma::colvec momentum_;
  const double momentum_coef_;
  const std::unique_ptr<Rcpp::Function> multiple_epochs_function_;
  arma::colvec objective_function_values_;
  arma::colvec objective_function_values_candidates_;
  double* objective_function_values_candidates_ptr_;
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
  const bool r_progress_;
  const unsigned int regression_response_count_;
  arma::mat result_coefficients_;
  arma::mat result_residuals_;
  double result_value_;
  std::unique_ptr<RProgress::RProgress> rProgress_;
  arma::mat segment_coefficients_;
  const int segment_count_;
  arma::colvec segment_indices_;
  unsigned int t = 1;
  const double trim_;
  const bool use_warm_start_;
  const double vanilla_percentage_;
  const arma::mat variance_estimate_;
  arma::mat warm_start_;
  friend fastcpd::test::FastcpdTest;
};

}  // namespace classes

}  // namespace fastcpd

#endif  // FASTCPD_CLASS_H_
