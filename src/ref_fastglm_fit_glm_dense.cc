// This is a modified copy of fastglm/src/fit_glm_dense.cpp

#define EIGEN_DONT_PARALLELIZE

#include "ref_fastglm_fit_glm_dense.h"

#include "ref_fastglm_fit_glm.h"
#include "ref_fastglm_glm.h"

using ::Eigen::MatrixXd;

typedef MatrixXd::Index Index;

List fit_glm(NumericMatrix Xs, NumericVector ys, NumericVector weightss,
             NumericVector offsets, NumericVector starts, NumericVector etas,
             int type, double tol, int maxit, string family) {
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
  NumericVector (*dev_resids)(const Map<VectorXd> &, const Eigen::VectorXd &,
                              const Map<VectorXd> &);
  if (family == "gaussian") {
    mu_eta = &mu_eta_gaussian;
    linkinv = &linkinv_gaussian;
    var = &var_gaussian;
    valideta = &valideta_gaussian;
    validmu = &validmu_gaussian;
    dev_resids = &dev_resids_gaussian;
  } else if (family == "binomial") {
    mu_eta = &mu_eta_binomial;
    linkinv = &linkinv_binomial;
    var = &var_binomial;
    valideta = &valideta_binomial;
    validmu = &validmu_binomial;
    dev_resids = &dev_resids_binomial;
  } else if (family == "poisson") {
    mu_eta = &mu_eta_poisson;
    linkinv = &linkinv_poisson;
    var = &var_poisson;
    valideta = &valideta_poisson;
    validmu = &validmu_poisson;
    dev_resids = &dev_resids_poisson;
  } else {
    throw invalid_argument("invalid family");
  }
  NumericVector mus = linkinv(as<Map<VectorXd>>(etas));
  const Map<VectorXd> mu_init(as<Map<VectorXd>>(mus));
  const Map<VectorXd> eta_init(as<Map<VectorXd>>(etas));
  Index n = X.rows();
  if ((Index)y.size() != n) throw invalid_argument("size mismatch");

  // instantiate fitting class
  GlmBase<Eigen::VectorXd, Eigen::MatrixXd> *glm_solver = NULL;

  bool is_big_matrix = false;

  glm_solver = new glm(X, y, weights, offset, var, mu_eta, linkinv, dev_resids,
                       valideta, validmu, tol, maxit, type, is_big_matrix);

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
  if (family == "binomial") {
    if ((mu.array() > 1 - eps).any() || (mu.array() < eps).any())
      warning("fit_glm: fitted probabilities numerically 0 or 1 occurred");
  }
  if (family == "poisson") {
    if ((mu.array() < eps).any())
      warning("fit_glm: fitted rates numerically 0 occurred");
  }

  return List::create(
      _["coefficients"] = beta, _["se"] = se, _["fitted.values"] = mu,
      _["linear.predictors"] = eta, _["deviance"] = dev, _["weights"] = wts,
      _["prior.weights"] = pweights, _["rank"] = rank, _["df.residual"] = df,
      _["residuals"] = as<NumericVector>(wrap(y - mu)) / mu_eta(eta),
      _["iter"] = iters, _["converged"] = converged);
}
