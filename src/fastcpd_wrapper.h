#ifndef FASTCPD_WRAPPER_H_
#define FASTCPD_WRAPPER_H_

#include "fastcpd_types.h"

namespace fastcpd::classes {

struct ColMat {
  mat data;

  operator colvec() const {
    // TODO(doccstat): Add a warning if the matrix has more than one column.
    return data.as_col();
  }

  operator mat() const {
    return data;
  }

  operator rowvec() const {
    // TODO(doccstat): Add a warning if the matrix has more than one column.
    return data.as_col().t();
  }
};

struct CostResult {
  ColMat par;
  ColMat residuals;
  double value;
};

}  // namespace fastcpd::classes

#endif  // FASTCPD_WRAPPER_H_
