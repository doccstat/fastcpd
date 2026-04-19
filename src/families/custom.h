#ifndef FASTCPD_FAMILIES_CUSTOM_H_
#define FASTCPD_FAMILIES_CUSTOM_H_

#include "fastcpd_family.h"

namespace fastcpd {
namespace families {

class CustomFamily : public BaseFamily {
 public:
  static constexpr bool has_optimized_run = false;
  static constexpr bool is_custom = true;

  template <typename T>
  static void GetNllPeltValue(T* instance, unsigned int const segment_start,
                              unsigned int const segment_end, bool const cv,
                              std::optional<arma::colvec> const& start) {
    GetNllPelt(instance, segment_start, segment_end, cv, start);
  }

  template <typename T>
  static void GetNllPelt(T* instance, unsigned int const segment_start,
                         unsigned int const segment_end, bool const cv,
                         std::optional<arma::colvec> const& start) {
    if (instance->cost_gradient_ || instance->cost_hessian_) {
      instance->GetOptimizedCostResult(segment_start, segment_end);
    } else {
      instance->result_coefficients_ = arma::colvec();
      instance->result_residuals_ = arma::mat();
      instance->result_value_ = instance->cost_function_pelt_(instance->data_.rows(segment_start, segment_end));
    }
  }

  template <typename T>
  static double GetNllSen(T* instance, unsigned int const segment_start,
                          unsigned int const segment_end,
                          arma::colvec const& theta) {
    return instance->cost_function_sen_(instance->data_.rows(segment_start, segment_end), theta);
  }

  template <typename T>
  static arma::colvec GetGradient(T* instance, unsigned int const segment_start,
                                  unsigned int const segment_end,
                                  arma::colvec const& theta) {
    return instance->cost_gradient_(instance->data_.rows(segment_start, segment_end), theta);
  }

  template <typename T>
  static arma::mat GetHessian(T* instance, unsigned int const segment_start,
                              unsigned int const segment_end,
                              arma::colvec const& theta) {
    return instance->cost_hessian_(instance->data_.rows(segment_start, segment_end), theta);
  }
};

} // namespace families
} // namespace fastcpd

#endif
