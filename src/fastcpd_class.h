#ifndef FASTCPD_CLASS_H_
#define FASTCPD_CLASS_H_

#include <memory>
#include <unordered_map>

// TODO(doccstat): Fix the order issue.
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
    arma::mat (Fastcpd::*hessian_)(const unsigned int segment_start,
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

  // Stop the clock and create an R object with `name`.
  void CreateRClock(const std::string name);

  // Initialize \code{theta_hat}, \code{theta_sum}, and \code{hessian}.
  void CreateSenParameters();

  // Initialize theta_hat_t_t to be the estimate in the segment.
  void CreateSegmentStatistics();

  // Initialize \code{theta_hat}, \code{theta_sum}, and \code{hessian} and
  // theta_hat_t_t to be the estimate in the segment.
  void CreateSegmentStatisticsAndSenParameters();

  // Set \code{theta_sum} for a specific column.
  void CreateThetaSum(const unsigned int col, arma::colvec new_theta_sum);

  double GetCostAdjustmentValue(const unsigned nrows);

  // Solve logistic/poisson regression using Gradient Descent Extension to the
  // multivariate case
  //
  // @param data A data frame containing the data to be segmented.
  // @param theta Estimate of the parameters. If null, the function will
  //   estimate the parameters.
  // @param lambda Lambda for L1 regularization. Only used for lasso.
  // @param cv Whether to perform cross-validation to find the best lambda.
  // @param start Starting point for the optimization for warm start.
  //   only used in mean change and lm.
  //
  // @return Negative log likelihood of the corresponding data with the given
  //   family.
  CostResult GetCostResult(const unsigned int segment_start,
                           const unsigned int segment_end,
                           Rcpp::Nullable<arma::colvec> theta,
                           const bool cv = false,
                           Rcpp::Nullable<arma::colvec> start = R_NilValue);
  Rcpp::List GetChangePointSet(const arma::colvec raw_cp_set);
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
  arma::colvec GetObjectiveFunctionValues(const arma::colvec& fvec,
                                          const arma::ucolvec& r_t_set,
                                          unsigned int r_t_count,
                                          unsigned int t);

  // Update \code{theta_hat}, \code{theta_sum}, and \code{hessian}.
  //
  // @param data_segment A data frame containing a segment of the data.
  //
  // @return A list containing new values of \code{theta_hat}, \code{theta_sum},
  //   and \code{hessian}.
  CostResult GetOptimizedCostResult(const unsigned int segment_start,
                                    const unsigned int segment_end);

  void Step(unsigned int t, arma::ucolvec& r_t_set, unsigned int& r_t_count,
            arma::colvec& cp_sets, arma::colvec& fvec);

  // Adjust cost value for MBIC and MDL.
  double UpdateCostValue(double value, const unsigned int nrows);

  arma::colvec UpdateChangePointSet(const arma::colvec raw_cp_set);

  // Append new values to \code{fastcpd_parameters}.
  void UpdateSenParameters(const unsigned int t);

  void UpdateSenParameters(const unsigned int segment_start,
                           const unsigned int segment_end,
                           const unsigned int i);

  void UpdateSenParametersStep(const int segment_start, const int segment_end,
                               const int i);

  // Update the cost values for the segmentation.
  //
  // @param data A data frame containing the data to be segmented.
  // @param tau Start of the current segment.
  // @param i Index of the current data in the whole data set.
  // @param family Family of the model.
  // @param momentum Momentum from the previous iteration.
  // @param epsilon Epsilon to avoid numerical issues. Only used for binomial
  //   and poisson.
  // @param line_search A vector containing the line search coefficients.
  //
  // @return A list containing new values of \code{theta_hat}, \code{theta_sum},
  //   \code{hessian}, and \code{momentum}.
  Rcpp::List UpdateSenParametersSteps(const int segment_start,
                                      const unsigned int segment_end,
                                      const int i, arma::colvec momentum);

  // Append a new slice to \code{hessian}.
  void UpdateHessian(arma::mat new_hessian);

  // Prune the slices of \code{hessian}.
  void UpdateHessian(arma::ucolvec pruned_left);

  // Update \code{hessian} for a specific slice.
  void UpdateHessian(const unsigned int slice, arma::mat new_hessian);

  // Update \code{momentum}.
  void UpdateMomentum(arma::colvec new_momentum);

  // Start the clock tick for `name`.
  void UpdateRClockTick(const std::string name);

  // Stop the clock tick for `name`.
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

  // `act_num_` is used in Lasso and Gaussian families only.
  arma::colvec act_num_;

  // `beta` is the initial cost value.
  double beta;

  // `cost` is the cost function to be used.
  const std::unique_ptr<Rcpp::Function> cost;

  // Adjustment to the cost function.
  const std::string cost_adjustment;

  // `cost_gradient` is the gradient of the cost function to be used.
  const std::unique_ptr<Rcpp::Function> cost_gradient;

  // `cost_hessian` is the Hessian of the cost function to be used.
  const std::unique_ptr<Rcpp::Function> cost_hessian;

  const bool cp_only_;

  // Dimension of the data.
  const unsigned int d;

  // `data` is the data set to be segmented.
  const arma::mat data;

  // The number of data points.
  const unsigned int data_n_rows_;

  // The number of data columns.
  const unsigned int data_n_cols_;

  // `epsilon` is the epsilon to avoid numerical issues. Only used for binomial
  // and poisson.
  const double epsilon;

  // `error_sd` is used in Gaussian family only.
  arma::colvec err_sd_;

  // `family` is the family of the model.
  const std::string family;

  // Rcpp::Function to calculate the gradient at the current data.
  //
  // @param data A data frame containing the data to be segmented.
  // @param theta Estimated theta from the previous iteration.
  // @param family Family of the model.
  // @param order Order of the time series models.
  //
  // @return Gradient at the current data.
  arma::colvec (Fastcpd::*get_gradient_)(const unsigned int segment_start,
                                         const unsigned int segment_end,
                                         const arma::colvec& theta);

  // Rcpp::Function to calculate the Hessian matrix at the current data.
  //
  // @param data A data frame containing the data to be segmented.
  // @param theta Estimated theta from the previous iteration.
  //
  // @return Hessian at the current data.
  arma::mat (Fastcpd::*get_hessian_)(const unsigned int segment_start,
                                     const unsigned int segment_end,
                                     const arma::colvec& theta);

  CostResult (Fastcpd::*get_nll_pelt_)(
      const unsigned int segment_start, const unsigned int segment_end,
      const bool cv, const Rcpp::Nullable<arma::colvec>& start);

  double (Fastcpd::*get_nll_sen_)(const unsigned int segment_start,
                                  const unsigned int segment_end,
                                  arma::colvec theta);

  // `hessian_` stores the Hessian matrix up to the current data point.
  arma::cube hessian_;

  // `lambda` is the lambda for L1 regularization. Only used for lasso. This
  // parameter stores the penalty without the segment length scaling.
  double lambda;
  arma::colvec line_search;

  // Lower bound of the parameters to be estimated during the optimization.
  const arma::colvec lower;

  // `min_idx` is the index of the minimum objective value.
  // This value is stored to avoid reallocation of memory.
  unsigned int min_idx;

  // `min_obj` is the minimum objective value.
  // This value is stored to avoid reallocation of memory.
  double min_obj;

  // Momentum will be used in the update step if `momentum_coef_` is not 0.
  arma::colvec momentum_;
  const double momentum_coef_;
  const std::unique_ptr<Rcpp::Function> multiple_epochs_function;
  const arma::colvec order;

  // `p` is the number of parameters to be estimated.
  const unsigned int p;

  // Number of response variables in regression.
  const unsigned int p_response;

  // Indices after pruning.
  arma::ucolvec pruned_left;
  unsigned int pruned_left_n_elem;
  const double pruning_coef;
  const std::string r_clock;
  const bool r_progress;
  Rcpp::Clock rClock;
  std::unique_ptr<RProgress::RProgress> rProgress;
  const int segment_count;
  arma::colvec segment_indices_;

  // Create a matrix to store the estimated coefficients in each segment,
  // where each row represents estimated coefficients for a segment.
  arma::mat segment_theta_hat_;

  // Matrix storing the warm starts.
  arma::mat start;

  // `theta_hat_` stores the estimated coefficients up to the current data
  // point.
  arma::mat theta_hat_;

  // `theta_sum` stores the sum of estimated coefficients up to the current data
  // point.
  arma::mat theta_sum_;
  const double trim;

  // Static mapping from family names to their function sets
  static const std::unordered_map<std::string, FunctionSet> family_function_map;

  // Upper bound of the parameters to be estimated during the optimization.
  const arma::colvec upper;
  const double vanilla_percentage;
  const arma::mat variance_estimate;
  const bool use_warm_start_;
  arma::mat zero_data_;
  unsigned int zero_data_n_cols_;
  unsigned int zero_data_n_rows_;
  double* zero_data_ptr_;
  friend fastcpd::test::FastcpdTest;
};

}  // namespace fastcpd::classes

#endif  // FASTCPD_CLASS_H_
