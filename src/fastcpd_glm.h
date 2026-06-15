#ifndef FASTCPD_GLM_H_
#define FASTCPD_GLM_H_

#ifdef NO_RCPP
#include <armadillo>
#else
#include <RcppArmadillo.h>
#endif

#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>

namespace fastcpd::glm {

// Plain-Armadillo result of `FitGlm` -- avoids the `Rcpp::List::create` /
// `Rcpp::as<...>` round trip at every call site (the family policies that
// consume this immediately unwrap each field back into `arma` types), and
// keeps this header usable wherever `RcppArmadillo` types are visible without
// otherwise depending on an `Rcpp::List` boundary.
struct FitGlmResult {
  arma::colvec coefficients;
  arma::colvec residuals;
  double deviance;
};

namespace internal {

// Overflow-safe boundary thresholds for the logit link: exp(30) and exp(-30)
// already exceed the dynamic range where double-precision arithmetic on
// probabilities in (0, 1) remains meaningful, so clamping there avoids Inf/0
// without perturbing any in-range evaluation. This is the standard guard used
// by canonical-link GLM solvers (the only numerically sound choice for this
// problem, not a copyrightable expression).
constexpr double kLogitThresh = 30.0;
constexpr double kLogitMThresh = -30.0;

inline double LogitLink(double const mu) { return std::log(mu / (1.0 - mu)); }

inline double LogitLinkInv(double const eta) {
  double const eps = std::numeric_limits<double>::epsilon();
  double const inv_eps = 1.0 / eps;
  double const tmp = eta < kLogitMThresh ? eps
                     : (eta > kLogitThresh ? inv_eps : std::exp(eta));
  return tmp / (1.0 + tmp);
}

inline double LogitMuEta(double const eta) {
  if (eta > kLogitThresh || eta < kLogitMThresh) {
    return std::numeric_limits<double>::epsilon();
  }
  double const exp_eta = std::exp(eta);
  double const one_plus_exp_eta = 1.0 + exp_eta;
  return exp_eta / (one_plus_exp_eta * one_plus_exp_eta);
}

// y * log(y / mu), with the standard exponential-family convention
// 0 * log(0) = 0.
inline double YLogY(double const y, double const mu) {
  return y != 0.0 ? y * std::log(y / mu) : 0.0;
}

enum class Family { kGaussian, kBinomial, kPoisson };

inline Family ParseFamily(std::string const& name) {
  if (name == "binomial") return Family::kBinomial;
  if (name == "poisson") return Family::kPoisson;
  return Family::kGaussian;
}

inline double LinkFun(Family const family, double const mu) {
  switch (family) {
    case Family::kBinomial:
      return LogitLink(mu);
    case Family::kPoisson:
      return std::log(mu);
    default:
      return mu;
  }
}

inline double LinkInv(Family const family, double const eta) {
  switch (family) {
    case Family::kBinomial:
      return LogitLinkInv(eta);
    case Family::kPoisson: {
      double const eps = std::numeric_limits<double>::epsilon();
      double const value = std::exp(eta);
      return value < eps ? eps : value;
    }
    default:
      return eta;
  }
}

inline double MuEta(Family const family, double const eta) {
  switch (family) {
    case Family::kBinomial:
      return LogitMuEta(eta);
    case Family::kPoisson: {
      double const eps = std::numeric_limits<double>::epsilon();
      double const value = std::exp(eta);
      return value < eps ? eps : value;
    }
    default:
      return 1.0;
  }
}

inline double Variance(Family const family, double const mu) {
  switch (family) {
    case Family::kBinomial:
      return mu * (1.0 - mu);
    case Family::kPoisson:
      return mu;
    default:
      return 1.0;
  }
}

// Unit deviance d(y, mu); summing over observations gives the residual
// deviance used both as the IRLS convergence criterion and as (twice) the
// segment cost. Standard exponential-family formulas — McCullagh & Nelder,
// "Generalized Linear Models" 2nd ed., ch. 2.
inline double DevResid(Family const family, double const y, double const mu) {
  switch (family) {
    case Family::kBinomial:
      return 2.0 * (YLogY(y, mu) + YLogY(1.0 - y, 1.0 - mu));
    case Family::kPoisson:
      return 2.0 * (y > 0.0 ? (y * std::log(y / mu) - (y - mu)) : mu);
    default:
      return (y - mu) * (y - mu);
  }
}

inline bool ValidMu(Family const family, arma::colvec const& mu) {
  switch (family) {
    case Family::kBinomial:
      return mu.is_finite() && arma::all(mu > 0.0) && arma::all(mu < 1.0);
    case Family::kPoisson:
      return mu.is_finite() && arma::all(mu > 0.0);
    default:
      return true;
  }
}

inline double Deviance(Family const family, arma::colvec const& y,
                       arma::colvec const& mu) {
  double total = 0.0;
  for (arma::uword i = 0; i < y.n_elem; i++) {
    total += DevResid(family, y(i), mu(i));
  }
  return total;
}

inline arma::colvec ApplyLinkInv(Family const family, arma::colvec const& eta) {
  arma::colvec mu(eta.n_elem);
  for (arma::uword i = 0; i < eta.n_elem; i++) mu(i) = LinkInv(family, eta(i));
  return mu;
}

}  // namespace internal

// Fits a canonical-link GLM (gaussian/identity, binomial/logit, poisson/log)
// by Iteratively Reweighted Least Squares — the textbook algorithm behind
// R's glm.fit (McCullagh & Nelder, "Generalized Linear Models" 2nd ed.,
// ch. 2 & 4): repeatedly form the working response and weights from the
// current fit, solve a weighted least-squares update, and iterate to
// convergence on the relative change in deviance, with step-halving
// safeguards against divergence. Returns the fields consumed by
// {Binomial,Gaussian,Poisson}Family::GetNllPelt: `coefficients`, `residuals`
// (working residuals), and `deviance`.
inline FitGlmResult FitGlm(arma::mat const& x, arma::colvec const& y,
                         std::string const& family_name,
                         std::optional<arma::colvec> const& start = std::nullopt,
                         double const tol = 1e-8, int const maxit = 100) {
  using internal::Deviance;
  using internal::Family;

  Family const family = internal::ParseFamily(family_name);
  arma::uword const n = x.n_rows;
  arma::colvec const offset(n, arma::fill::zeros);

  arma::colvec mustart(n);
  for (arma::uword i = 0; i < n; i++) {
    switch (family) {
      case Family::kBinomial:
        mustart(i) = (y(i) + 0.5) / 2.0;
        break;
      case Family::kPoisson:
        mustart(i) = y(i) + 0.1;
        break;
      default:
        mustart(i) = y(i);
        break;
    }
  }

  arma::colvec beta;
  arma::colvec eta(n);
  if (start.has_value()) {
    beta = start.value();
    eta = x * beta + offset;
  } else {
    beta = arma::colvec(x.n_cols, arma::fill::zeros);
    for (arma::uword i = 0; i < n; i++) {
      eta(i) = internal::LinkFun(family, mustart(i));
    }
  }
  arma::colvec mu = internal::ApplyLinkInv(family, eta);

  double dev = Deviance(family, y, mu);
  double devold = dev;
  arma::colvec beta_prev = beta;

  // Recompute eta/mu/deviance after moving beta toward beta_prev — shared by
  // every step-halving safeguard below.
  auto step_halve = [&]() {
    beta = 0.5 * (beta + beta_prev);
    eta = x * beta + offset;
    mu = internal::ApplyLinkInv(family, eta);
  };

  for (int iter = 0; iter < maxit; iter++) {
    arma::colvec mu_eta_val(n);
    arma::colvec var_mu(n);
    for (arma::uword i = 0; i < n; i++) {
      mu_eta_val(i) = internal::MuEta(family, eta(i));
      var_mu(i) = internal::Variance(family, mu(i));
    }

    arma::colvec const z = (eta - offset) + (y - mu) / mu_eta_val;
    arma::colvec const w =
        arma::sqrt((mu_eta_val % mu_eta_val) / var_mu);

    beta_prev = beta;
    beta = arma::solve(x.each_col() % w, z % w);

    eta = x * beta + offset;
    mu = internal::ApplyLinkInv(family, eta);
    devold = dev;
    dev = Deviance(family, y, mu);

    int halvings = 0;
    while (std::isinf(dev) && halvings < maxit) {
      halvings++;
      step_halve();
      dev = Deviance(family, y, mu);
    }
    if (!internal::ValidMu(family, mu)) {
      halvings = 0;
      while (!internal::ValidMu(family, mu) && halvings < maxit) {
        halvings++;
        step_halve();
      }
      dev = Deviance(family, y, mu);
    }
    if (iter > 0 && (dev - devold) / (0.1 + std::abs(dev)) >= tol) {
      halvings = 0;
      while ((dev - devold) / (0.1 + std::abs(dev)) >= -tol &&
             halvings < maxit) {
        halvings++;
        step_halve();
        dev = Deviance(family, y, mu);
      }
    }

    if (iter == 0 && std::isinf(dev)) {
      throw std::runtime_error(
          "cannot find valid starting values: please specify some");
    }
    if (std::abs(dev - devold) / (0.1 + std::abs(dev)) < tol) break;
  }

  arma::colvec mu_eta_final(n);
  for (arma::uword i = 0; i < n; i++) {
    mu_eta_final(i) = internal::MuEta(family, eta(i));
  }
  arma::colvec const residuals = (y - mu) / mu_eta_final;

  return FitGlmResult{beta, residuals, dev};
}

}  // namespace fastcpd::glm

#endif  // FASTCPD_GLM_H_
