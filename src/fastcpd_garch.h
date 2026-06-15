#ifndef FASTCPD_GARCH_H_
#define FASTCPD_GARCH_H_

#ifdef NO_RCPP
#include <armadillo>
#else
#include <RcppArmadillo.h>
#endif

#include <algorithm>
#include <cmath>
#include <utility>

#include "fastcpd_optim.h"

namespace fastcpd::garch {

// Plain-Armadillo result of `FitGarch` -- avoids the `Rcpp::List::create` /
// `Rcpp::as<...>` round trip at the call site (`GarchFamily::GetNllPelt`
// immediately unwraps each field back into `arma`/`double`), keeping this
// header's only Rcpp dependency the `RcppArmadillo` include itself.
struct FitGarchResult {
  arma::colvec coef;
  double n_likeli;
  arma::colvec residuals;
};

namespace internal {

// Sentinel objective value for parameter vectors outside the admissible
// region (a_0 > 0, a_i, b_j >= 0 — the standard non-negativity constraint
// that keeps every conditional variance positive). Mirrors the common
// "return a large value" handling of constraints inside unconstrained
// quasi-Newton solvers.
constexpr double kInvalidPenalty = 1.0e10;

// Tolerance for the non-negativity constraint check below. A boundary-
// constrained optimum (e.g. an ARCH(q) model nested inside GARCH(p, q) with
// some b_j == 0) is approached by the unconstrained BFGS iteration from
// both sides; the line search can overshoot the boundary by an amount on
// the order of machine epsilon (e.g. -1e-16). Treating such points as
// "invalid" would return the penalty value *and* a zero gradient, which
// makes the gradient-norm stopping criterion fire immediately at a garbage,
// penalized point. Accepting points within kBoundaryTolerance of the
// boundary keeps the objective and gradient smooth across it.
constexpr double kBoundaryTolerance = 1.0e-8;

// Negative log-likelihood and analytical gradient of a GARCH(p, q) model
// (Bollerslev 1986, "Generalized Autoregressive Conditional
// Heteroskedasticity", Journal of Econometrics 31). Apart from an additive
// constant and the presample contribution:
//
//   h_t  = a_0 + sum_{j=1}^q a_j x_{t-j}^2 + sum_{j=1}^p b_j h_{t-j}
//   NLL  = (1/2) sum_{t >= max(p,q)} [log(h_t) + x_t^2 / h_t]
//
// The gradient follows from the recursive chain rule on h_t:
//   dh_t/da_0 = 1         + sum_{j=1}^p b_j dh_{t-j}/da_0
//   dh_t/da_k = x_{t-k}^2 + sum_{j=1}^p b_j dh_{t-j}/da_k
//   dh_t/db_k = h_{t-k}   + sum_{j=1}^p b_j dh_{t-j}/db_k
//   dNLL/dtheta_k = sum_t [(1/2)(1 - x_t^2/h_t)/h_t] * dh_t/dtheta_k
//
// Presample values h_t, t < max(p, q), are backcast to the empirical
// uncentred second moment mean(x^2) — the standard fixed presample used when
// maximising the GARCH likelihood (Bollerslev 1986, sec. 3-4); since it does
// not depend on theta, its derivative contribution is zero.
inline std::pair<double, arma::colvec> NegLogLikelihood(
    arma::colvec const& x, int const p, int const q,
    arma::colvec const& theta) {
  int const n = static_cast<int>(x.n_elem);
  int const ncoef = static_cast<int>(theta.n_elem);
  int const max_pq = std::max(p, q);

  if (theta(0) <= 0.0 ||
      arma::any(theta.tail(ncoef - 1) < -kBoundaryTolerance)) {
    return {kInvalidPenalty, arma::colvec(ncoef, arma::fill::zeros)};
  }

  double const presample_h = arma::dot(x, x) / n;

  arma::colvec h(n, arma::fill::zeros);
  arma::mat dh(n, ncoef, arma::fill::zeros);
  for (int t = 0; t < max_pq && t < n; t++) {
    h(t) = presample_h;
    dh(t, 0) = 1.0;
  }

  double nll = 0.0;
  arma::colvec grad(ncoef, arma::fill::zeros);

  for (int t = max_pq; t < n; t++) {
    double h_t = theta(0);
    for (int j = 1; j <= q; j++) h_t += theta(j) * x(t - j) * x(t - j);
    for (int j = 1; j <= p; j++) h_t += theta(q + j) * h(t - j);
    h(t) = h_t;

    double dh_a0 = 1.0;
    for (int j = 1; j <= p; j++) dh_a0 += theta(q + j) * dh(t - j, 0);
    dh(t, 0) = dh_a0;

    for (int k = 1; k <= q; k++) {
      double dh_ak = x(t - k) * x(t - k);
      for (int j = 1; j <= p; j++) dh_ak += theta(q + j) * dh(t - j, k);
      dh(t, k) = dh_ak;
    }
    for (int k = 1; k <= p; k++) {
      double dh_bk = h(t - k);
      for (int j = 1; j <= p; j++) dh_bk += theta(q + j) * dh(t - j, q + k);
      dh(t, q + k) = dh_bk;
    }

    double const x2 = x(t) * x(t);
    nll += std::log(h_t) + x2 / h_t;

    double const dl_dh = 0.5 * (1.0 - x2 / h_t) / h_t;
    for (int k = 0; k < ncoef; k++) grad(k) += dl_dh * dh(t, k);
  }

  return {0.5 * nll, grad};
}

// Conditional-variance recursion used for prediction/residuals. The presample
// is backcast to the model's own stationary variance
// omega / (1 - sum(alpha) - sum(beta)) — the textbook unconditional variance
// of a covariance-stationary GARCH(p, q) process (Bollerslev 1986, eq. 6) —
// rather than the empirical second moment used while fitting, matching the
// "non-genuine" in-sample prediction convention.
inline arma::colvec PredictConditionalVariance(arma::colvec const& x,
                                               int const p, int const q,
                                               arma::colvec const& theta) {
  int const n = static_cast<int>(x.n_elem);
  int const max_pq = std::max(p, q);
  double const sum_coef = arma::accu(theta.tail(theta.n_elem - 1));
  double const stationary_var = theta(0) / (1.0 - sum_coef);

  arma::colvec h(n);
  for (int t = 0; t < max_pq && t < n; t++) h(t) = stationary_var;
  for (int t = max_pq; t < n; t++) {
    double h_t = theta(0);
    for (int j = 1; j <= q; j++) h_t += theta(j) * x(t - j) * x(t - j);
    for (int j = 1; j <= p; j++) h_t += theta(q + j) * h(t - j);
    h(t) = h_t;
  }
  return h;
}

}  // namespace internal

// Fits a GARCH(p, q) model by maximum likelihood (Bollerslev 1986) using a
// from-scratch BFGS optimiser driven by the analytical gradient — see
// internal::NegLogLikelihood and fastcpd::optim::BfgsMinimize — replacing the
// vendored PORT `dsumsl`-based fit. Returns the fields consumed by
// GarchFamily::GetNllPelt: `coef` (fitted parameters), `n_likeli` (negative
// log-likelihood at the optimum), and `residuals` (standardised residuals
// x_t / sqrt(h_t), NaN for the first max(p, q) presample points).
inline FitGarchResult FitGarch(arma::colvec const& x, arma::colvec const& order) {
  int const p = static_cast<int>(order(0));
  int const q = static_cast<int>(order(1));
  int const ncoef = p + q + 1;
  int const n = static_cast<int>(x.n_elem);
  int const max_pq = std::max(p, q);

  // Starting values follow the convention used by standard GARCH fitting
  // routines: small positive ARCH/GARCH coefficients, with the intercept
  // chosen so the implied stationary variance matches the unconditional
  // variance of x. For a (G)ARCH process E[x] = 0, so the unconditional
  // variance equals the uncentred second moment E[x^2] — estimated here by
  // mean(x^2), matching the `presample_h` backcast used inside
  // NegLogLikelihood. Unlike the centred sample variance arma::var(x), this
  // estimator is well-defined (finite, non-negative) for segments as short
  // as a single observation, which PELT routinely evaluates as candidates
  // enter the pruned set; arma::var(x) divides by (n - 1) and produces NaN
  // for n == 1, poisoning theta(0) and the subsequent fit.
  constexpr double kSmall = 0.05;
  double const mean_sq = arma::dot(x, x) / static_cast<double>(x.n_elem);
  arma::colvec theta(ncoef, arma::fill::value(kSmall));
  theta(0) = mean_sq * (1.0 - kSmall * (ncoef - 1));

  // NegLogLikelihood sums n per-observation terms, so both its value and its
  // gradient scale with n — an *absolute* gradient-norm stopping criterion
  // (as used by BfgsMinimize) would then need a different grad_tol for every
  // segment length. Minimising the mean per-observation NLL instead (an
  // equivalent problem: scaling the objective by the positive constant 1/n
  // does not change the argmin) keeps the gradient at order O(1) regardless
  // of n, so a single fixed grad_tol is meaningful for every candidate
  // segment PELT evaluates.
  double const n_obs = static_cast<double>(n);
  auto objective = [&](arma::colvec const& params) {
    auto [f, grad] = internal::NegLogLikelihood(x, p, q, params);
    return std::pair<double, arma::colvec>(f / n_obs, grad / n_obs);
  };
  theta = fastcpd::optim::BfgsMinimize(
      objective, theta, arma::zeros<arma::colvec>(ncoef),
      arma::colvec(ncoef).fill(arma::datum::inf));

  // Project the converged point onto the closed feasible region
  // {a_0 > 0, a_i, b_j >= 0}: BFGS approaches boundary-constrained optima
  // (e.g. a vanishing GARCH coefficient) from an unconstrained iteration and
  // can overshoot by machine-epsilon (see kBoundaryTolerance above). Clamping
  // here both reports clean coefficients (exact 0 rather than -1e-16) and
  // ensures the final objective/residual evaluation below sees a feasible
  // point.
  theta = arma::clamp(theta, 0.0, arma::datum::inf);

  double const nll = internal::NegLogLikelihood(x, p, q, theta).first;
  arma::colvec const h = internal::PredictConditionalVariance(x, p, q, theta);

  arma::colvec residuals(n);
  for (int t = 0; t < n; t++) {
    residuals(t) = (t < max_pq) ? arma::datum::nan : x(t) / std::sqrt(h(t));
  }

  return FitGarchResult{theta, nll, residuals};
}

}  // namespace fastcpd::garch

#endif  // FASTCPD_GARCH_H_
