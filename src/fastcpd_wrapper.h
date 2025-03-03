#ifndef FASTCPD_WRAPPER_H_
#define FASTCPD_WRAPPER_H_

// fastcpd_wrapper.h
//
// This header file defines wrapper types for fastcpd.
//
// It provides conversion operators for matrix types using RcppArmadillo.

#include "RcppArmadillo.h"

namespace fastcpd {
namespace classes {

// A lightweight wrapper around an Armadillo matrix.
// Provides conversion operators to common Armadillo vector and matrix types.
struct ColMat {
  arma::mat data;

  // Converts to an Armadillo column vector.
  // TODO(doccstat): Add a warning if the matrix has more than one column.
  operator arma::colvec() const { return data.as_col(); }

  // Converts to an Armadillo matrix.
  operator arma::mat() const { return data; }

  // Converts to an Armadillo row vector.
  // TODO(doccstat): Add a warning if the matrix has more than one column.
  operator arma::rowvec() const { return data.as_col().t(); }
};

// A structure to hold the result of a cost calculation.
struct CostResult {
  ColMat par;
  ColMat residuals;
  double value;
};

}  // namespace classes
}  // namespace fastcpd

#endif  // FASTCPD_WRAPPER_H_
