// This is a modified copy of fastglm/inst/include/glm_base.h

#ifndef GLM_BASE_H
#define GLM_BASE_H

#include <RcppArmadillo.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace std;

using Eigen::Map;

template <typename VecTypeX, typename MatTypeX>
class GlmBase {
 protected:
  int const nvars;  // dimension of beta
  int const nobs;   // number of rows

  VecTypeX beta;       // parameters to be optimized
  VecTypeX beta_prev;  // auxiliary parameters

  VecTypeX eta;
  VecTypeX var_mu;
  VecTypeX mu_eta;
  VecTypeX mu;
  VecTypeX z;
  VecTypeX w;
  MatTypeX vcov;
  VecTypeX se;
  double dev, devold, devnull;

  int maxit;   // max iterations
  double tol;  // tolerance for convergence
  bool conv;

  virtual bool converged() {
    if (std::abs(dev - devold) / (0.1 + std::abs(dev)) < tol) {
      return true;
    } else {
      return false;
    }
  }

  virtual void update_eta() {}

  virtual void update_var_mu() {}

  virtual void update_mu_eta() {}

  virtual void update_mu() {}

  virtual void update_z() {}

  virtual void update_w() {}

  virtual void step_halve() {}

  virtual void run_step_halving(int& iterr) {}

  virtual void update_dev_resids() {}

  virtual void update_dev_resids_dont_update_old() {}

  virtual void solve_wls(int iter) {}

  virtual void save_se() {}

 public:
  GlmBase(int n_, int p_, double tol_ = 1e-6, int maxit_ = 100)
      : nvars(p_),
        nobs(n_),
        beta(p_),
        beta_prev(p_),  // allocate space but do not set values
        eta(n_),
        var_mu(n_),
        mu_eta(n_),
        mu(n_),
        z(n_),
        w(n_),
        vcov(p_, p_),
        se(p_),
        maxit(maxit_),
        tol(tol_) {}

  virtual ~GlmBase() {}

  virtual void init_parms(Map<VecTypeX> const& start_, Map<VecTypeX> const& mu_,
                          Map<VecTypeX> const& eta_) {}

  void update_beta() {
    // VecTypeX newbeta(nvars);
    next_beta(beta);
    // beta.swap(newbeta);
  }

  int solve(int maxit) {
    int i;

    conv = false;

    for (i = 0; i < maxit; ++i) {
      update_var_mu();

      update_mu_eta();

      update_z();

      update_w();

      solve_wls(i);

      update_eta();

      update_mu();

      update_dev_resids();

      run_step_halving(i);

      if (std::isinf(dev) && i == 0) {
        stop("cannot find valid starting values: please specify some");
      }

      if (converged()) {
        conv = true;
        break;
      }
    }

    save_se();

    return std::min(i + 1, maxit);
  }

  virtual VecTypeX get_beta() { return beta; }
  virtual VecTypeX get_eta() { return eta; }
  virtual VecTypeX get_se() { return se; }
  virtual VecTypeX get_mu() { return mu; }
  virtual VecTypeX get_weights() { return w; }
  virtual VecTypeX get_w() { return w.array().square(); }
  virtual double get_dev() { return dev; }
  virtual int get_rank() { return nvars; }
  virtual MatTypeX get_vcov() { return vcov; }
  virtual bool get_converged() { return conv; }
};

#endif  // GLM_BASE_H
