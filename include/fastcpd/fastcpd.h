#ifndef FASTCPD_FASTCPD_H_
#define FASTCPD_FASTCPD_H_

#include <armadillo>

#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <string>
#include <vector>

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

struct VarianceArmaRow {
  unsigned int order;
  double sigma2;
  double aic;
  double bic;
};

struct VarianceArmaResult {
  std::vector<VarianceArmaRow> table;
  double sigma2_aic;
  double sigma2_bic;
};

struct ConfidenceInterval {
  std::string parm;
  std::string method;
  std::string bootstrap;
  unsigned int index = 0;
  unsigned int segment = 0;
  unsigned int parameter = 0;
  double estimate = std::numeric_limits<double>::quiet_NaN();
  double lower = std::numeric_limits<double>::quiet_NaN();
  double upper = std::numeric_limits<double>::quiet_NaN();
  double detection_rate = std::numeric_limits<double>::quiet_NaN();
  double profile_min = std::numeric_limits<double>::quiet_NaN();
  double cutoff = std::numeric_limits<double>::quiet_NaN();
  double se = std::numeric_limits<double>::quiet_NaN();
  double level = 0.95;
};

struct ConfidenceOptions {
  std::string parm = "cp";
  std::string method;
  double level = 0.95;
  unsigned int bootstrap_replicates = 999;
  std::string bootstrap = "nonparametric";
  std::optional<unsigned int> window;
  unsigned int min_segment_length = 2;
  std::optional<std::int32_t> seed;
  // C++ results intentionally do not copy the input or full fit options.
  // Supply the options used for fitting when refit/profile behavior depends
  // on non-default detector controls.
  Options detector_options;
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

arma::mat estimate_variance_mean(arma::mat const& data);
double estimate_variance_median(arma::mat const& data);
arma::mat estimate_variance_linear_regression(
    arma::mat const& data, unsigned int d = 1, unsigned int block_size = 0,
    double outlier_iqr = std::numeric_limits<double>::infinity());
arma::mat estimate_variance_lm(
    arma::mat const& data, unsigned int d = 1, unsigned int block_size = 0,
    double outlier_iqr = std::numeric_limits<double>::infinity());
VarianceArmaResult estimate_variance_arma(arma::colvec const& data,
                                          unsigned int p, unsigned int q,
                                          unsigned int max_order = 0);

std::vector<ConfidenceInterval> confint(
    Result const& result, arma::mat const& data,
    ConfidenceOptions options = {});
std::vector<ConfidenceInterval> confint(
    Result const& result, arma::colvec const& data,
    ConfidenceOptions options = {});

}  // namespace fastcpd

#endif  // FASTCPD_FASTCPD_H_
