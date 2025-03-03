#ifndef FASTCPD_TEST_H_
#define FASTCPD_TEST_H_

#include "fastcpd_wrapper.h"

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
                          const unsigned int segment_end, arma::colvec theta);

  // Computes the negative log-likelihood for the PELT model.
  static fastcpd::classes::CostResult GetNllPelt(
      const arma::mat& data, const unsigned int segment_start,
      const unsigned int segment_end, const bool cv,
      const Rcpp::Nullable<arma::colvec>& start);

  // Updates the theta sum.
  static arma::mat UpdateThetaSum(const unsigned int col,
                                  arma::colvec old_theta_sum,
                                  arma::colvec new_theta_sum);
};

}  // namespace test
}  // namespace fastcpd

#endif  // FASTCPD_TEST_H_
