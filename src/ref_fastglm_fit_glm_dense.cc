// This is a modified copy of fastglm/src/fit_glm_dense.cpp

#define EIGEN_DONT_PARALLELIZE

#include <Rcpp.h>
#include "ref_fastglm_glm.h"
#include <RcppEigen.h>

using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef MatrixXd::Index Index;

// extern SEXP logit_mu_eta(SEXP eta);
extern "C"
{
  SEXP logit_link(SEXP mu);
  SEXP logit_linkinv(SEXP eta);
  SEXP logit_mu_eta(SEXP eta);
  SEXP binomial_dev_resids(SEXP y, SEXP mu, SEXP wt);
}

bool valideta_gaussian(const Eigen::VectorXd &eta)
{
  return true;
}

bool valideta_binomial(const Eigen::VectorXd &eta)
{
  return true;
}

bool valideta_poisson(const Eigen::VectorXd &eta)
{
  return true;
}

bool validmu_gaussian(const Eigen::VectorXd &mu)
{
  return true;
}

bool validmu_binomial(const Eigen::VectorXd &mu)
{
  return mu.allFinite() && (mu.array() > 0).all() && (mu.array() < 1).all();
}

bool validmu_poisson(const Eigen::VectorXd &mu)
{
  return mu.allFinite() && (mu.array() > 0).all();
}

// Gaussian deviance residuals: wt * ((y - mu)^2)
Rcpp::NumericVector dev_resids_gaussian(const Map<VectorXd> &y, const VectorXd &mu, const Map<VectorXd> &wt)
{
  int n = y.size();
  Rcpp::NumericVector ans(n);
  for (int i = 0; i < n; i++)
  {
    ans[i] = wt[i] * std::pow(y[i] - mu[i], 2);
  }
  return ans;
}

// Binomial deviance residuals: simply call the exported C function.
// Note: binomial_dev_resids_cpp is assumed to wrap the SEXP version defined elsewhere.
Rcpp::NumericVector dev_resids_binomial(const Map<VectorXd> &y, const VectorXd &mu, const Map<VectorXd> &wt)
{
  // Convert the Eigen vectors to Rcpp NumericVectors.
  Rcpp::NumericVector R_y = Rcpp::wrap(y);
  Rcpp::NumericVector R_mu = Rcpp::wrap(mu);
  Rcpp::NumericVector R_wt = Rcpp::wrap(wt);

  // Call the exported function that returns SEXP.
  SEXP res = binomial_dev_resids(R_y, R_mu, R_wt);
  return Rcpp::NumericVector(res);
}

// Poisson deviance residuals:
//   r = mu * wt,
//   for indices where y > 0, set r = wt * (y * log(y/mu) - (y - mu))
//   and return 2 * r.
Rcpp::NumericVector dev_resids_poisson(const Map<VectorXd> &y, const VectorXd &mu, const Map<VectorXd> &wt)
{
  int n = y.size();
  Rcpp::NumericVector ans(n);
  for (int i = 0; i < n; i++)
  {
    double r = mu[i] * wt[i];
    if (y[i] > 0)
    {
      r = wt[i] * (y[i] * std::log(y[i] / mu[i]) - (y[i] - mu[i]));
    }
    ans[i] = 2.0 * r;
  }
  return ans;
}

// Gaussian variance: always 1 for each element.
NumericVector var_gaussian(const VectorXd &mu)
{
  int n = mu.size();
  NumericVector ans(n);
  std::fill(ans.begin(), ans.end(), 1.0);
  return ans;
}

// Binomial variance: mu * (1 - mu)
NumericVector var_binomial(const VectorXd &mu)
{
  int n = mu.size();
  NumericVector ans(n);
  for (int i = 0; i < n; i++)
  {
    ans[i] = mu[i] * (1.0 - mu[i]);
  }
  return ans;
}

// Poisson variance: just mu
NumericVector var_poisson(const VectorXd &mu)
{
  int n = mu.size();
  NumericVector ans(n);
  for (int i = 0; i < n; i++)
  {
    ans[i] = mu[i];
  }
  return ans;
}

// Gaussian link inverse: identity (eta).
NumericVector linkinv_gaussian(const VectorXd &eta)
{
  // return NumericVector(1);
  // Simply wrap the Eigen vector as a NumericVector.
  return Rcpp::wrap(eta);
}

// Binomial link inverse: delegate to logit_linkinv_cpp.
NumericVector linkinv_binomial(const VectorXd &eta)
{
  // return NumericVector(1);
  // Wrap eta into a NumericVector.
  NumericVector R_eta = Rcpp::wrap(eta);
  // Call the exported function.
  SEXP res = logit_linkinv(R_eta);
  return NumericVector(res);
}

// Poisson link inverse: pmax(exp(eta), .Machine$double.eps)
NumericVector linkinv_poisson(const VectorXd &eta)
{
  // return NumericVector(1);
  int n = eta.size();
  NumericVector ans(n);
  // .Machine$double.eps in R is typically 2.220446e-16.
  // We'll use the C++ equivalent:
  double eps = std::numeric_limits<double>::epsilon();
  for (int i = 0; i < n; i++)
  {
    double value = std::exp(eta[i]);
    ans[i] = (value < eps) ? eps : value;
  }
  return ans;
}

// Gaussian mu.eta: returns a vector of ones, analogous to rep.int(1, length(eta))
NumericVector mu_eta_gaussian(const VectorXd &eta)
{
  return NumericVector(eta.size(), 1.0);
}

// Binomial mu.eta: delegate to the exported logit_mu_eta_cpp function.
// It is assumed that logit_mu_eta_cpp is declared elsewhere and available.
NumericVector mu_eta_binomial(const VectorXd &eta)
{
  NumericVector R_eta = wrap(eta);
  SEXP res = logit_mu_eta(R_eta);
  return NumericVector(res);
}

// Poisson mu.eta: computes pmax(exp(eta), .Machine$double.eps)
NumericVector mu_eta_poisson(const VectorXd &eta)
{
  int n = eta.size();
  NumericVector ans(n);
  // Use the C++ equivalent of .Machine$double.eps.
  double eps = std::numeric_limits<double>::epsilon();
  for (int i = 0; i < n; i++)
  {
    double value = std::exp(eta[i]);
    ans[i] = (value < eps) ? eps : value;
  }
  return ans;
}

NumericVector linkfun_gaussian(const NumericVector &mu)
{
  return mu;
}

NumericVector linkfun_binomial(const NumericVector &mu)
{
  SEXP res = logit_link(wrap(mu));
  return as<NumericVector>(res);
}

NumericVector linkfun_poisson(const NumericVector &mu)
{
  return Rcpp::log(mu);
}

List fit_glm(Rcpp::NumericMatrix Xs,
             Rcpp::NumericVector ys,
             Rcpp::NumericVector weightss,
             Rcpp::NumericVector offsets,
             Rcpp::NumericVector starts,
             Rcpp::NumericVector etas,
             int type,
             double tol,
             int maxit,
             std::string family)
{
  const Map<MatrixXd> X(as<Map<MatrixXd>>(Xs));
  const Map<VectorXd> y(as<Map<VectorXd>>(ys));
  const Map<VectorXd> weights(as<Map<VectorXd>>(weightss));
  const Map<VectorXd> offset(as<Map<VectorXd>>(offsets));
  const Map<VectorXd> beta_init(as<Map<VectorXd>>(starts));
  NumericVector (*mu_eta)(const Eigen::VectorXd &);
  NumericVector (*linkinv)(const Eigen::VectorXd &);
  NumericVector (*var)(const Eigen::VectorXd &);
  bool (*valideta)(const Eigen::VectorXd &);
  bool (*validmu)(const Eigen::VectorXd &);
  Rcpp::NumericVector (*dev_resids)(const Map<VectorXd> &, const Eigen::VectorXd &, const Map<VectorXd> &);
  if (family == "gaussian")
  {
    mu_eta = &mu_eta_gaussian;
    linkinv = &linkinv_gaussian;
    var = &var_gaussian;
    valideta = &valideta_gaussian;
    validmu = &validmu_gaussian;
    dev_resids = &dev_resids_gaussian;
  }
  else if (family == "binomial")
  {
    mu_eta = &mu_eta_binomial;
    linkinv = &linkinv_binomial;
    var = &var_binomial;
    valideta = &valideta_binomial;
    validmu = &validmu_binomial;
    dev_resids = &dev_resids_binomial;
  }
  else if (family == "poisson")
  {
    mu_eta = &mu_eta_poisson;
    linkinv = &linkinv_poisson;
    var = &var_poisson;
    valideta = &valideta_poisson;
    validmu = &validmu_poisson;
    dev_resids = &dev_resids_poisson;
  }
  else
  {
    throw invalid_argument("invalid family");
  }
  NumericVector mus = linkinv(as<Map<VectorXd>>(etas));
  const Map<VectorXd> mu_init(as<Map<VectorXd>>(mus));
  const Map<VectorXd> eta_init(as<Map<VectorXd>>(etas));
  Index n = X.rows();
  if ((Index)y.size() != n)
    throw invalid_argument("size mismatch");

  // instantiate fitting class
  GlmBase<Eigen::VectorXd, Eigen::MatrixXd> *glm_solver = NULL;

  bool is_big_matrix = false;

  glm_solver = new glm(X, y, weights, offset,
                       var, mu_eta, linkinv, dev_resids,
                       valideta, validmu, tol, maxit, type,
                       is_big_matrix);

  // initialize parameters
  glm_solver->init_parms(beta_init, mu_init, eta_init);

  // maximize likelihood
  int iters = glm_solver->solve(maxit);

  VectorXd beta = glm_solver->get_beta();
  VectorXd se = glm_solver->get_se();
  VectorXd mu = glm_solver->get_mu();
  VectorXd eta = glm_solver->get_eta();
  VectorXd wts = glm_solver->get_w();
  VectorXd pweights = glm_solver->get_weights();

  double dev = glm_solver->get_dev();
  int rank = glm_solver->get_rank();
  bool converged = glm_solver->get_converged();

  int df = X.rows() - rank;

  delete glm_solver;

  double eps = 10 * std::numeric_limits<double>::epsilon();
  if (family == "binomial")
  {
    if ((mu.array() > 1 - eps).any() || (mu.array() < eps).any())
      warning("fit_glm: fitted probabilities numerically 0 or 1 occurred");
  }
  if (family == "poisson")
  {
    if ((mu.array() < eps).any())
      warning("fit_glm: fitted rates numerically 0 occurred");
  }

  return List::create(_["coefficients"] = beta,
                      _["se"] = se,
                      _["fitted.values"] = mu,
                      _["linear.predictors"] = eta,
                      _["deviance"] = dev,
                      _["weights"] = wts,
                      _["prior.weights"] = pweights,
                      _["rank"] = rank,
                      _["df.residual"] = df,
                      _["residuals"] = as<NumericVector>(wrap(y - mu)) / mu_eta(eta),
                      _["iter"] = iters,
                      _["converged"] = converged);
}

Eigen::MatrixXd colMax_dense(const Eigen::Map<Eigen::MatrixXd> &A)
{
  Eigen::VectorXd colM = A.colwise().maxCoeff();
  return colM;
}

Eigen::MatrixXd colMin_dense(const Eigen::Map<Eigen::MatrixXd> &A)
{
  Eigen::VectorXd colM = A.colwise().minCoeff();
  return colM;
}

List fastglm(NumericMatrix x,
             SEXP y,
             std::string family,
             Nullable<NumericVector> start,
             Nullable<NumericVector> weights,
             Nullable<NumericVector> offset,
             Nullable<NumericVector> etastart,
             Nullable<NumericVector> mustart,
             int method,
             double tol,
             int maxit)
{
  // Determine number of observations and whether y is a matrix.
  int nobs;
  bool yIsMatrix = false;
  NumericMatrix yMat = as<NumericMatrix>(y);
  NumericVector yVec;
  if (yMat.ncol() > 1)
  {
    nobs = yMat.nrow();
    yIsMatrix = true;
  }
  else
  {
    yVec = yMat(_, 0);
    nobs = yVec.size();
  }

  // Set default weights (rep(1, nobs)) and offset (rep(0, nobs)) if not provided.
  NumericVector wt = weights.isNotNull() ? as<NumericVector>(weights)
                                         : NumericVector(nobs, 1.0);
  NumericVector off = offset.isNotNull() ? as<NumericVector>(offset)
                                         : NumericVector(nobs, 0.0);

  // Save a copy of mustart if provided.
  Nullable<NumericVector> mukeep = mustart;

  // Variables to hold computed "mustart" and an auxiliary vector n.
  NumericVector must;
  NumericVector n;

  if (family == "binomial")
  {
    if (!yIsMatrix)
    {
      // Single-column y: assume numeric (ignoring factor conversion).
      n = NumericVector(nobs, 1.0);
      // For observations with zero weight, set y to 0.
      for (int i = 0; i < nobs; i++)
      {
        if (wt[i] == 0)
          yVec[i] = 0;
      }
      // Check that y is in [0,1].
      for (int i = 0; i < nobs; i++)
      {
        if (yVec[i] < 0 || yVec[i] > 1)
          stop("y values must be 0 <= y <= 1");
      }
      must = NumericVector(nobs);
      for (int i = 0; i < nobs; i++)
      {
        must[i] = (wt[i] * yVec[i] + 0.5) / (wt[i] + 1);
      }
      // Warn if weighted successes are not integer.
      for (int i = 0; i < nobs; i++)
      {
        double m_val = wt[i] * yVec[i];
        if (std::abs(m_val - std::round(m_val)) > 0.001)
          warning("non-integer #successes in a binomial glm!");
      }
    }
    else
    {
      // y is a two-column matrix.
      if (yMat.ncol() != 2)
        stop("For binomial family, y must have 1 or 2 columns");
      // Warn if counts are non-integer.
      for (int i = 0; i < yMat.nrow(); i++)
      {
        for (int j = 0; j < yMat.ncol(); j++)
        {
          if (std::abs(yMat(i, j) - std::round(yMat(i, j))) > 0.001)
            warning("non-integer counts in a binomial glm!");
        }
      }
      NumericVector y1 = yMat(_, 0); // first column
      NumericVector y2 = yMat(_, 1); // second column
      n = NumericVector(nobs);
      for (int i = 0; i < nobs; i++)
      {
        n[i] = y1[i] + y2[i];
      }
      yVec = NumericVector(nobs);
      for (int i = 0; i < nobs; i++)
      {
        if (n[i] == 0)
          yVec[i] = 0;
        else
          yVec[i] = y1[i] / n[i];
      }
      // Adjust weights.
      for (int i = 0; i < nobs; i++)
      {
        wt[i] = wt[i] * n[i];
      }
      must = NumericVector(nobs);
      for (int i = 0; i < nobs; i++)
      {
        must[i] = (n[i] * yVec[i] + 0.5) / (n[i] + 1);
      }
    }
  }
  else if (family == "poisson")
  {
    if (yIsMatrix)
      stop("Poisson family expects y to be a vector");
    for (int i = 0; i < nobs; i++)
    {
      if (yVec[i] < 0)
        stop("negative values not allowed for the 'Poisson' family");
    }
    n = NumericVector(nobs, 1.0);
    must = NumericVector(nobs);
    for (int i = 0; i < nobs; i++)
    {
      must[i] = yVec[i] + 0.1;
    }
  }
  else if (family == "gaussian")
  {
    if (yIsMatrix)
      stop("Gaussian family expects y to be a vector");
    n = NumericVector(nobs, 1.0);
    must = NumericVector(nobs);
    for (int i = 0; i < nobs; i++)
    {
      must[i] = yVec[i];
    }
  }
  else
  {
    stop("Unsupported family");
  }

  // If mustart was provided, override the computed value.
  if (mukeep.isNotNull())
    must = as<NumericVector>(mukeep);

  // Define the link function as a lambda.
  // For gaussian, the identity; for binomial, call logit_link_cpp; for poisson, log.
  NumericVector (*linkfun)(const NumericVector &);
  if (family == "gaussian")
  {
    linkfun = &linkfun_gaussian;
  }
  else if (family == "binomial")
  {
    linkfun = &linkfun_binomial;
  }
  else if (family == "poisson")
  {
    linkfun = &linkfun_poisson;
  }
  else
  {
    stop("Unsupported family");
  }

  // Ensure y is numeric: if y was a matrix, we now use the processed yVec.
  // Otherwise, yVec is already set.

  // Compute eta: use etastart if provided; else, if start is provided compute offset + x %*% start;
  // otherwise, use linkfun(must).
  NumericVector eta;
  if (etastart.isNotNull())
  {
    eta = as<NumericVector>(etastart);
  }
  else if (start.isNotNull())
  {
    NumericVector startVec = as<NumericVector>(start);
    if (startVec.size() != x.ncol())
      stop("Length of start must equal number of columns of x");
    int nrow = x.nrow();
    NumericVector xstart(nrow, 0.0);
    for (int i = 0; i < nrow; i++)
    {
      double sum = 0.0;
      for (int j = 0; j < x.ncol(); j++)
      {
        sum += x(i, j) * startVec[j];
      }
      xstart[i] = sum;
    }
    eta = off; // start with offset
    for (int i = 0; i < nrow; i++)
    {
      eta[i] += xstart[i];
    }
  }
  else
  {
    eta = linkfun(must);
  }

  // If start is not provided, default to a zero vector of appropriate length.
  NumericVector startVec;
  if (start.isNotNull())
  {
    startVec = as<NumericVector>(start);
  }
  else
  {
    startVec = NumericVector(x.ncol(), 0.0);
  }

  // Finally, call fastglm with the processed arguments.
  // Note: we use yVec as the final response.
  List res = fit_glm(x, yVec, wt, off, startVec, eta, method, tol, maxit, family);
  return res;
}
