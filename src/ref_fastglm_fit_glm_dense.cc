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
