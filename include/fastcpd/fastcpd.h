#ifndef FASTCPD_FASTCPD_H_
#define FASTCPD_FASTCPD_H_

#include <armadillo>

#include <cstdint>
#include <functional>
#include <optional>
#include <string>

namespace fastcpd {

using CostFunctionPelt = std::function<double(arma::mat const&)>;
using CostFunctionSen =
    std::function<double(arma::mat const&, arma::colvec const&)>;
using CostGradientFunction =
    std::function<arma::colvec(arma::mat const&, arma::colvec const&)>;
using CostHessianFunction =
    std::function<arma::mat(arma::mat const&, arma::colvec const&)>;
using MultipleEpochsFunction = std::function<unsigned int(unsigned int)>;

struct Result {
  arma::colvec raw_change_points;
  arma::colvec change_points;
  arma::colvec cost_values;
  arma::mat residuals;
  arma::mat thetas;
  std::string family;
  arma::colvec order;
  bool cp_only = false;
};

struct Options {
  std::string family = "mean";
  std::optional<double> beta;
  std::string beta_criterion = "MBIC";
  // Empty selects the wrapper default: MBIC generally and BIC for KCP.
  std::string cost_adjustment;
  bool cp_only = false;
  double epsilon = 1e-10;
  arma::colvec line_search = arma::colvec{1.0};
  arma::colvec lower;
  arma::colvec upper;
  double momentum_coef = 0.0;
  MultipleEpochsFunction multiple_epochs =
      [](unsigned int) -> unsigned int { return 0u; };
  arma::colvec order = arma::colvec{0.0, 0.0, 0.0};
  int p = 0;
  unsigned int p_response = 0;
  std::optional<double> pruning_coef;
  int segment_count = 10;
  double trim = 0.0;
  double vanilla_percentage = 0.0;
  arma::mat variance_estimate;
  bool warm_start = false;
  bool show_progress = false;
  bool include_mean = false;
  std::optional<std::int32_t> seed;

  CostFunctionPelt cost_pelt;
  CostFunctionSen cost_sen;
  CostGradientFunction cost_gradient;
  CostHessianFunction cost_hessian;
};

Result detect(arma::mat const& data, Options options = {});
Result detect(arma::colvec const& data, Options options = {});

Result detect_mean(arma::mat const& data, Options options = {});
Result detect_variance(arma::mat const& data, Options options = {});
Result detect_meanvariance(arma::mat const& data, Options options = {});
Result detect_mean_variance(arma::mat const& data, Options options = {});
Result detect_exponential(arma::mat const& data, Options options = {});
Result detect_lm(arma::mat const& data, Options options = {});
Result detect_linear_regression(arma::mat const& data, Options options = {});
Result detect_lasso(arma::mat const& data, Options options = {});
Result detect_binomial(arma::mat const& data, Options options = {});
Result detect_logistic_regression(arma::mat const& data,
                                  Options options = {});
Result detect_poisson(arma::mat const& data, Options options = {});
Result detect_poisson_regression(arma::mat const& data,
                                 Options options = {});
Result detect_quantile(arma::mat const& data, Options options = {});
Result detect_quantile_regression(arma::mat const& data,
                                  Options options = {});
Result detect_ar(arma::colvec const& data, Options options = {});
Result detect_arma(arma::colvec const& data, Options options = {});
Result detect_arima(arma::colvec const& data, Options options = {});
Result detect_garch(arma::colvec const& data, Options options = {});
Result detect_var(arma::mat const& data, Options options = {});
Result detect_mgaussian(arma::mat const& data, Options options = {});
Result detect_rank(arma::mat const& data, Options options = {});
Result detect_kernel(arma::mat const& data, Options options = {});
Result detect_kcp(arma::mat const& data, Options options = {});

Result mean(arma::mat const& data, Options options = {});
Result variance(arma::mat const& data, Options options = {});
Result meanvariance(arma::mat const& data, Options options = {});
Result exponential(arma::mat const& data, Options options = {});
Result gaussian(arma::mat const& data, Options options = {});
Result lm(arma::mat const& data, Options options = {});
Result lasso(arma::mat const& data, Options options = {});
Result binomial(arma::mat const& data, Options options = {});
Result poisson(arma::mat const& data, Options options = {});
Result quantile(arma::mat const& data, Options options = {});
Result ar(arma::colvec const& data, Options options = {});
Result arima(arma::colvec const& data, Options options = {});
Result var(arma::mat const& data, Options options = {});
Result mgaussian(arma::mat const& data, Options options = {});
Result rank(arma::mat const& data, Options options = {});
Result kernel(arma::mat const& data, Options options = {});
Result kcp(arma::mat const& data, Options options = {});

}  // namespace fastcpd

#endif  // FASTCPD_FASTCPD_H_
