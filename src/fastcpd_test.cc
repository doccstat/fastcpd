// fastcpd_test.cc
//
// Implementation of the FastcpdTest helper class for testing fastcpd
// functionality.

#include "fastcpd.h"

namespace fastcpd {
namespace test {

//------------------------------------------------------------------------------
// GetGradientArma
// Computes the gradient using the Armadillo implementation.
//------------------------------------------------------------------------------
arma::colvec FastcpdTest::GetGradientArma(arma::mat const& data,
                                          unsigned int const segment_start,
                                          unsigned int const segment_end,
                                          arma::colvec const& theta) {
  fastcpd::classes::Fastcpd fastcpd_instance(
      /* beta */ 0,
      /* cost */ std::nullopt,
      /* cost_pelt */ nullptr,
      /* cost_sen */ nullptr,
      /* cost_adjustment */ "MBIC",
      /* cost_gradient */ nullptr,
      /* cost_hessian */ nullptr,
      /* cp_only */ true,
      /* data */ data,
      /* epsilon */ 0.0,
      /* family */ "arma",
      /* multiple_epochs_function */ nullptr,
      /* line_search */ arma::colvec(),
      /* lower */ arma::colvec(),
      /* momentum_coef */ 0.0,
      /* order */ arma::colvec({3, 2}),
      /* p */ 0,
      /* p_response */ 0,
      /* pruning_coef */ 0,
      /* r_progress */ false,
      /* segment_count */ 0,
      /* trim */ 0,
      /* upper */ arma::colvec(),
      /* vanilla_percentage */ 0.0,
      /* variance_estimate */ arma::mat(),
      /* warm_start */ false);
  return fastcpd_instance.GetGradientArma(segment_start, segment_end, theta);
}

//------------------------------------------------------------------------------
// GetHessianArma
// Computes the Hessian using the Armadillo implementation.
//------------------------------------------------------------------------------
arma::mat FastcpdTest::GetHessianArma(arma::mat const& data,
                                      unsigned int const segment_start,
                                      unsigned int const segment_end,
                                      arma::colvec const& theta) {
  fastcpd::classes::Fastcpd fastcpd_instance(
      /* beta */ 0,
      /* cost */ std::nullopt,
      /* cost_pelt */ nullptr,
      /* cost_sen */ nullptr,
      /* cost_adjustment */ "MBIC",
      /* cost_gradient */ nullptr,
      /* cost_hessian */ nullptr,
      /* cp_only */ true,
      /* data */ data,
      /* epsilon */ 0.0,
      /* family */ "arma",
      /* multiple_epochs_function */ nullptr,
      /* line_search */ arma::colvec(),
      /* lower */ arma::colvec(),
      /* momentum_coef */ 0.0,
      /* order */ arma::colvec({3, 2}),
      /* p */ 0,
      /* p_response */ 0,
      /* pruning_coef */ 0,
      /* r_progress */ false,
      /* segment_count */ 0,
      /* trim */ 0,
      /* upper */ arma::colvec(),
      /* vanilla_percentage */ 0.0,
      /* variance_estimate */ arma::mat(),
      /* warm_start */ false);
  return fastcpd_instance.GetHessianArma(segment_start, segment_end, theta);
}

//------------------------------------------------------------------------------
// GetHessianBinomial
// Computes the Hessian for a binomial model.
//------------------------------------------------------------------------------
arma::mat FastcpdTest::GetHessianBinomial(arma::mat const& data,
                                          unsigned int const segment_start,
                                          unsigned int const segment_end,
                                          arma::colvec const& theta) {
  fastcpd::classes::Fastcpd fastcpd_instance(
      /* beta */ 0,
      /* cost */ std::nullopt,
      /* cost_pelt */ nullptr,
      /* cost_sen */ nullptr,
      /* cost_adjustment */ "MBIC",
      /* cost_gradient */ nullptr,
      /* cost_hessian */ nullptr,
      /* cp_only */ true,
      /* data */ data,
      /* epsilon */ 0.0,
      /* family */ "binomial",
      /* multiple_epochs_function */ nullptr,
      /* line_search */ arma::colvec(),
      /* lower */ arma::colvec(),
      /* momentum_coef */ 0.0,
      /* order */ arma::colvec(),
      /* p */ 0,
      /* p_response */ 0,
      /* pruning_coef */ 0,
      /* r_progress */ false,
      /* segment_count */ 0,
      /* trim */ 0,
      /* upper */ arma::colvec(),
      /* vanilla_percentage */ 0.0,
      /* variance_estimate */ arma::mat(),
      /* warm_start */ false);
  return fastcpd_instance.GetHessianBinomial(segment_start, segment_end, theta);
}

//------------------------------------------------------------------------------
// GetHessianPoisson
// Computes the Hessian for a Poisson model.
//------------------------------------------------------------------------------
arma::mat FastcpdTest::GetHessianPoisson(arma::mat const& data,
                                         unsigned int const segment_start,
                                         unsigned int const segment_end,
                                         arma::colvec const& theta) {
  fastcpd::classes::Fastcpd fastcpd_instance(
      /* beta */ 0,
      /* cost */ std::nullopt,
      /* cost_pelt */ nullptr,
      /* cost_sen */ nullptr,
      /* cost_adjustment */ "MBIC",
      /* cost_gradient */ nullptr,
      /* cost_hessian */ nullptr,
      /* cp_only */ true,
      /* data */ data,
      /* epsilon */ 0.0,
      /* family */ "poisson",
      /* multiple_epochs_function */ nullptr,
      /* line_search */ arma::colvec(),
      /* lower */ arma::colvec(),
      /* momentum_coef */ 0.0,
      /* order */ arma::colvec(),
      /* p */ 0,
      /* p_response */ 0,
      /* pruning_coef */ 0,
      /* r_progress */ false,
      /* segment_count */ 0,
      /* trim */ 0,
      /* upper */ arma::colvec(),
      /* vanilla_percentage */ 0.0,
      /* variance_estimate */ arma::mat(),
      /* warm_start */ false);
  return fastcpd_instance.GetHessianPoisson(segment_start, segment_end, theta);
}

//------------------------------------------------------------------------------
// GetNllSen
// Computes the negative log-likelihood for the SEN model.
//------------------------------------------------------------------------------
double FastcpdTest::GetNllSen(arma::mat const& data,
                              unsigned int const segment_start,
                              unsigned int const segment_end,
                              arma::colvec const& theta) {
  fastcpd::classes::Fastcpd fastcpd_instance(
      /* beta */ 0,
      /* cost */ std::nullopt,
      /* cost_pelt */ nullptr,
      /* cost_sen */ nullptr,
      /* cost_adjustment */ "MBIC",
      /* cost_gradient */ nullptr,
      /* cost_hessian */ nullptr,
      /* cp_only */ true,
      /* data */ data,
      /* epsilon */ 0.0,
      /* family */ "arma",
      /* multiple_epochs_function */ nullptr,
      /* line_search */ arma::colvec(),
      /* lower */ arma::colvec(),
      /* momentum_coef */ 0.0,
      /* order */ arma::colvec({3, 2}),
      /* p */ 0,
      /* p_response */ 0,
      /* pruning_coef */ 0,
      /* r_progress */ false,
      /* segment_count */ 0,
      /* trim */ 0,
      /* upper */ arma::colvec(),
      /* vanilla_percentage */ 0.0,
      /* variance_estimate */ arma::mat(),
      /* warm_start */ false);
  return (fastcpd_instance.*fastcpd_instance.get_nll_sen_)(segment_start,
                                                           segment_end, theta);
}

//------------------------------------------------------------------------------
// GetNllPelt
// Computes the negative log-likelihood for the PELT model.
//------------------------------------------------------------------------------
std::tuple<arma::colvec, arma::mat, double> FastcpdTest::GetNllPelt(
    arma::mat const& data, unsigned int const segment_start,
    unsigned int const segment_end, bool const cv,
    std::optional<arma::colvec> const& start) {
  fastcpd::classes::Fastcpd fastcpd_instance(
      /* beta */ 0,
      /* cost */ std::nullopt,
      /* cost_pelt */ nullptr,
      /* cost_sen */ nullptr,
      /* cost_adjustment */ "MBIC",
      /* cost_gradient */ nullptr,
      /* cost_hessian */ nullptr,
      /* cp_only */ true,
      /* data */ data,
      /* epsilon */ 0.0,
      /* family */ "arma",
      /* multiple_epochs_function */ nullptr,
      /* line_search */ arma::colvec(),
      /* lower */ arma::colvec(),
      /* momentum_coef */ 0.0,
      /* order */ arma::colvec({3, 2}),
      /* p */ 0,
      /* p_response */ 0,
      /* pruning_coef */ 0,
      /* r_progress */ false,
      /* segment_count */ 0,
      /* trim */ 0,
      /* upper */ arma::colvec(),
      /* vanilla_percentage */ 0.0,
      /* variance_estimate */ arma::mat(),
      /* warm_start */ false);
  (fastcpd_instance.*fastcpd_instance.get_nll_pelt_)(segment_start, segment_end,
                                                     cv, start);
  return std::make_tuple(fastcpd_instance.result_coefficients_,
                         fastcpd_instance.result_residuals_,
                         fastcpd_instance.result_value_);
}

}  // namespace test
}  // namespace fastcpd
