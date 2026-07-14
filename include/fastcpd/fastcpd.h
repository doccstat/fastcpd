#ifndef FASTCPD_FASTCPD_H_
#define FASTCPD_FASTCPD_H_

#include <armadillo>

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
};

struct Options {
  std::string family = "mean";
  std::optional<double> beta;
  std::string beta_criterion = "MBIC";
  std::string cost_adjustment = "MBIC";
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

  CostFunctionPelt cost_pelt;
  CostFunctionSen cost_sen;
  CostGradientFunction cost_gradient;
  CostHessianFunction cost_hessian;
};

Result detect(arma::mat const& data, Options options = {});
Result detect(arma::colvec const& data, Options options = {});

Result mean(arma::mat const& data, Options options = {});
Result variance(arma::mat const& data, Options options = {});
Result meanvariance(arma::mat const& data, Options options = {});
Result exponential(arma::mat const& data, Options options = {});
Result gaussian(arma::mat const& data, Options options = {});

}  // namespace fastcpd

#endif  // FASTCPD_FASTCPD_H_
