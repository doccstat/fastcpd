#ifndef WRAPPERS_H_
#define WRAPPERS_H_

#include "fastcpd_types.h"

namespace fastcpd::wrappers {

class CostFunction {
 public:
  CostFunction(Function cost);

  List operator()(
      mat data,
      Nullable<colvec> theta,
      string family,  // UNUSED
      double lambda,  // UNUSED
      bool cv,  // UNUSED
      Nullable<colvec> start  // UNUSED
  );

 private:
  Function cost;
};

class CostGradient {
 public:
  CostGradient(Function cost_gradient);

  colvec operator()(
      mat data,
      colvec theta,
      string family  // UNUSED
  );

 private:
  Function cost_gradient;
};

class CostHessian {
 public:
  CostHessian(Function cost_hessian);

  mat operator()(
      mat data,
      colvec theta,
      string family,  // UNUSED
      double min_prob  // UNUSED
  );

 private:
  Function cost_hessian;
};

}  // namespace fastcpd::wrappers

#endif  // WRAPPERS_H_
