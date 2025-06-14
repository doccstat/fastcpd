// This is a modified copy of fastglm/R/fit_glm.R

#include "ref_fastglm_fit_glm.h"

#include "ref_fastglm_fit_glm_dense.h"

using ::Eigen::Map;
using ::Eigen::VectorXd;
using ::Rcpp::_;
using ::Rcpp::as;
using ::Rcpp::NumericMatrix;
using ::Rcpp::NumericVector;
using ::Rcpp::stop;
using ::Rcpp::warning;
using ::Rcpp::wrap;
using ::std::abs;
using ::std::exp;
using ::std::fill;
using ::std::pow;
using ::std::round;
using ::std::string;

// extern SEXP logit_mu_eta(SEXP eta);
extern "C" {
SEXP logit_link(SEXP mu);
SEXP logit_linkinv(SEXP eta);
SEXP logit_mu_eta(SEXP eta);
SEXP binomial_dev_resids(SEXP y, SEXP mu, SEXP wt);
}

bool valideta_gaussian(const VectorXd &eta) { return true; }

bool valideta_binomial(const VectorXd &eta) { return true; }

bool valideta_poisson(const VectorXd &eta) { return true; }

bool validmu_gaussian(const VectorXd &mu) { return true; }

bool validmu_binomial(const VectorXd &mu) {
  return mu.allFinite() && (mu.array() > 0).all() && (mu.array() < 1).all();
}

bool validmu_poisson(const VectorXd &mu) {
  return mu.allFinite() && (mu.array() > 0).all();
}

NumericVector dev_resids_gaussian(const Map<VectorXd> &y, const VectorXd &mu,
                                  const Map<VectorXd> &wt) {
  int n = y.size();
  NumericVector ans(n);
  for (int i = 0; i < n; i++) {
    ans[i] = wt[i] * pow(y[i] - mu[i], 2);
  }
  return ans;
}

NumericVector dev_resids_binomial(const Map<VectorXd> &y, const VectorXd &mu,
                                  const Map<VectorXd> &wt) {
  // Convert the Eigen vectors to Rcpp NumericVectors.
  NumericVector R_y = wrap(y);
  NumericVector R_mu = wrap(mu);
  NumericVector R_wt = wrap(wt);

  // Call the exported function that returns SEXP.
  SEXP res = binomial_dev_resids(R_y, R_mu, R_wt);
  return NumericVector(res);
}

NumericVector dev_resids_poisson(const Map<VectorXd> &y, const VectorXd &mu,
                                 const Map<VectorXd> &wt) {
  int n = y.size();
  NumericVector ans(n);
  for (int i = 0; i < n; i++) {
    double r = mu[i] * wt[i];
    if (y[i] > 0) {
      r = wt[i] * (y[i] * std::log(y[i] / mu[i]) - (y[i] - mu[i]));
    }
    ans[i] = 2.0 * r;
  }
  return ans;
}

NumericVector var_gaussian(const VectorXd &mu) {
  int n = mu.size();
  NumericVector ans(n);
  fill(ans.begin(), ans.end(), 1.0);
  return ans;
}

NumericVector var_binomial(const VectorXd &mu) {
  int n = mu.size();
  NumericVector ans(n);
  for (int i = 0; i < n; i++) {
    ans[i] = mu[i] * (1.0 - mu[i]);
  }
  return ans;
}

NumericVector var_poisson(const VectorXd &mu) {
  int n = mu.size();
  NumericVector ans(n);
  for (int i = 0; i < n; i++) {
    ans[i] = mu[i];
  }
  return ans;
}

NumericVector linkinv_gaussian(const VectorXd &eta) {
  // return NumericVector(1);
  // Simply wrap the Eigen vector as a NumericVector.
  return wrap(eta);
}

NumericVector linkinv_binomial(const VectorXd &eta) {
  // return NumericVector(1);
  // Wrap eta into a NumericVector.
  NumericVector R_eta = wrap(eta);
  // Call the exported function.
  SEXP res = logit_linkinv(R_eta);
  return NumericVector(res);
}

NumericVector linkinv_poisson(const VectorXd &eta) {
  // return NumericVector(1);
  int n = eta.size();
  NumericVector ans(n);
  // .Machine$double.eps in R is typically 2.220446e-16.
  // We'll use the C++ equivalent:
  double eps = std::numeric_limits<double>::epsilon();
  for (int i = 0; i < n; i++) {
    double value = exp(eta[i]);
    ans[i] = (value < eps) ? eps : value;
  }
  return ans;
}

NumericVector mu_eta_gaussian(const VectorXd &eta) {
  return NumericVector(eta.size(), 1.0);
}

NumericVector mu_eta_binomial(const VectorXd &eta) {
  NumericVector R_eta = wrap(eta);
  SEXP res = logit_mu_eta(R_eta);
  return NumericVector(res);
}

NumericVector mu_eta_poisson(const VectorXd &eta) {
  int n = eta.size();
  NumericVector ans(n);
  // Use the C++ equivalent of .Machine$double.eps.
  double eps = std::numeric_limits<double>::epsilon();
  for (int i = 0; i < n; i++) {
    double value = exp(eta[i]);
    ans[i] = (value < eps) ? eps : value;
  }
  return ans;
}

NumericVector linkfun_gaussian(const NumericVector &mu) { return mu; }

NumericVector linkfun_binomial(const NumericVector &mu) {
  SEXP res = logit_link(wrap(mu));
  return as<NumericVector>(res);
}

NumericVector linkfun_poisson(const NumericVector &mu) { return Rcpp::log(mu); }

List fastglm(const mat &x_, const colvec &y, const string &family,
             std::optional<colvec> start, Nullable<NumericVector> weights,
             Nullable<NumericVector> offset, Nullable<NumericVector> etastart,
             Nullable<NumericVector> mustart, int method, double tol,
             int maxit) {
  NumericMatrix x = wrap(x_);
  // Determine number of observations and whether y is a matrix.
  int nobs;
  bool yIsMatrix = false;
  NumericMatrix yMat = wrap(y);
  NumericVector yVec;
  if (yMat.ncol() > 1) {
    nobs = yMat.nrow();
    yIsMatrix = true;
  } else {
    yVec = yMat(_, 0);
    nobs = yVec.size();
  }

  // Set default weights (rep(1, nobs)) and offset (rep(0, nobs)) if not
  // provided.
  NumericVector wt = weights.isNotNull() ? as<NumericVector>(weights)
                                         : NumericVector(nobs, 1.0);
  NumericVector off =
      offset.isNotNull() ? as<NumericVector>(offset) : NumericVector(nobs, 0.0);

  // Save a copy of mustart if provided.
  Nullable<NumericVector> mukeep = mustart;

  // Variables to hold computed "mustart" and an auxiliary vector n.
  NumericVector must;
  NumericVector n;

  if (family == "binomial") {
    if (!yIsMatrix) {
      // Single-column y: assume numeric (ignoring factor conversion).
      n = NumericVector(nobs, 1.0);
      // For observations with zero weight, set y to 0.
      for (int i = 0; i < nobs; i++) {
        if (wt[i] == 0) yVec[i] = 0;
      }
      // Check that y is in [0,1].
      for (int i = 0; i < nobs; i++) {
        if (yVec[i] < 0 || yVec[i] > 1) stop("y values must be 0 <= y <= 1");
      }
      must = NumericVector(nobs);
      for (int i = 0; i < nobs; i++) {
        must[i] = (wt[i] * yVec[i] + 0.5) / (wt[i] + 1);
      }
      // Warn if weighted successes are not integer.
      for (int i = 0; i < nobs; i++) {
        double m_val = wt[i] * yVec[i];
        if (abs(m_val - round(m_val)) > 0.001)
          warning("non-integer #successes in a binomial glm!");
      }
    } else {
      // y is a two-column matrix.
      if (yMat.ncol() != 2)
        stop("For binomial family, y must have 1 or 2 columns");
      // Warn if counts are non-integer.
      for (int i = 0; i < yMat.nrow(); i++) {
        for (int j = 0; j < yMat.ncol(); j++) {
          if (abs(yMat(i, j) - round(yMat(i, j))) > 0.001)
            warning("non-integer counts in a binomial glm!");
        }
      }
      NumericVector y1 = yMat(_, 0);  // first column
      NumericVector y2 = yMat(_, 1);  // second column
      n = NumericVector(nobs);
      for (int i = 0; i < nobs; i++) {
        n[i] = y1[i] + y2[i];
      }
      yVec = NumericVector(nobs);
      for (int i = 0; i < nobs; i++) {
        if (n[i] == 0)
          yVec[i] = 0;
        else
          yVec[i] = y1[i] / n[i];
      }
      // Adjust weights.
      for (int i = 0; i < nobs; i++) {
        wt[i] = wt[i] * n[i];
      }
      must = NumericVector(nobs);
      for (int i = 0; i < nobs; i++) {
        must[i] = (n[i] * yVec[i] + 0.5) / (n[i] + 1);
      }
    }
  } else if (family == "poisson") {
    if (yIsMatrix) stop("Poisson family expects y to be a vector");
    for (int i = 0; i < nobs; i++) {
      if (yVec[i] < 0)
        stop("negative values not allowed for the 'Poisson' family");
    }
    n = NumericVector(nobs, 1.0);
    must = NumericVector(nobs);
    for (int i = 0; i < nobs; i++) {
      must[i] = yVec[i] + 0.1;
    }
  } else if (family == "gaussian") {
    if (yIsMatrix) stop("Gaussian family expects y to be a vector");
    n = NumericVector(nobs, 1.0);
    must = NumericVector(nobs);
    for (int i = 0; i < nobs; i++) {
      must[i] = yVec[i];
    }
  } else {
    stop("Unsupported family");
  }

  // If mustart was provided, override the computed value.
  if (mukeep.isNotNull()) must = as<NumericVector>(mukeep);

  // Define the link function as a lambda.
  // For gaussian, the identity; for binomial, call logit_link_cpp; for poisson,
  // log.
  NumericVector (*linkfun)(const NumericVector &);
  if (family == "gaussian") {
    linkfun = &linkfun_gaussian;
  } else if (family == "binomial") {
    linkfun = &linkfun_binomial;
  } else if (family == "poisson") {
    linkfun = &linkfun_poisson;
  } else {
    stop("Unsupported family");
  }

  // Ensure y is numeric: if y was a matrix, we now use the processed yVec.
  // Otherwise, yVec is already set.

  // Compute eta: use etastart if provided; else, if start is provided compute
  // offset + x %*% start; otherwise, use linkfun(must).
  NumericVector eta;
  if (etastart.isNotNull()) {
    eta = as<NumericVector>(etastart);
  } else if (start.has_value()) {
    NumericVector startVec = wrap(start.value());
    if (startVec.size() != x.ncol())
      stop("Length of start must equal number of columns of x");
    int nrow = x.nrow();
    NumericVector xstart(nrow, 0.0);
    for (int i = 0; i < nrow; i++) {
      double sum = 0.0;
      for (int j = 0; j < x.ncol(); j++) {
        sum += x(i, j) * startVec[j];
      }
      xstart[i] = sum;
    }
    eta = off;  // start with offset
    for (int i = 0; i < nrow; i++) {
      eta[i] += xstart[i];
    }
  } else {
    eta = linkfun(must);
  }

  // If start is not provided, default to a zero vector of appropriate length.
  NumericVector startVec;
  if (start.has_value()) {
    startVec = wrap(start.value());
  } else {
    startVec = NumericVector(x.ncol(), 0.0);
  }

  // Finally, call fastglm with the processed arguments.
  // Note: we use yVec as the final response.
  List res =
      fit_glm(x, yVec, wt, off, startVec, eta, method, tol, maxit, family);
  return res;
}
