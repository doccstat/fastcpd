#ifndef FASTCPD_OPTIM_H_
#define FASTCPD_OPTIM_H_

#ifdef NO_RCPP
#include <armadillo>
#else
#include <RcppArmadillo.h>
#endif

#include <cmath>
#include <limits>

namespace fastcpd::optim {

// Brent's method for bound-constrained 1-D minimisation (Forsythe, Malcolm &
// Moler, "Computer Methods for Mathematical Computations", 1977, procedure
// `fmin` -- the algorithm behind R's `optimize()` / `optim(method = "Brent")`,
// reproduced here so the gradient/Hessian-driven PELT warm start
// (`Fastcpd::GetOptimizedCostResult`, `parameters_count_ == 1` branch) no
// longer needs to call back into R via `stats::optim`). `tol` matches
// `optim`'s default `reltol = sqrt(.Machine$double.eps)`.
template <typename Objective>
double BrentMinimize(Objective const& objective, double const lower,
                     double const upper,
                     double const tol = 1.4901161193847656e-08) {
  constexpr double kSquaredInverseGoldenRatio = 0.3819660112501051;
  double const eps = std::sqrt(std::numeric_limits<double>::epsilon());
  double const tol3 = tol / 3.0;

  double a = lower, b = upper;
  double v = a + kSquaredInverseGoldenRatio * (b - a);
  double w = v, x = v;
  double d = 0.0, e = 0.0;
  double fx = objective(x);
  double fv = fx, fw = fx;

  for (;;) {
    double const xm = (a + b) * 0.5;
    double const tol1 = eps * std::abs(x) + tol3;
    double const t2 = tol1 * 2.0;

    if (std::abs(x - xm) <= t2 - (b - a) * 0.5) break;

    double p = 0.0, q = 0.0, r = 0.0;
    if (std::abs(e) > tol1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2.0;
      if (q > 0.0) p = -p; else q = -q;
      r = e;
      e = d;
    }

    if (std::abs(p) >= std::abs(q * 0.5 * r) || p <= q * (a - x) ||
        p >= q * (b - x)) {
      e = (x < xm) ? (b - x) : (a - x);
      d = kSquaredInverseGoldenRatio * e;
    } else {
      d = p / q;
      double const u = x + d;
      if (u - a < t2 || b - u < t2) {
        d = (x >= xm) ? -tol1 : tol1;
      }
    }

    double const u =
        (std::abs(d) >= tol1) ? (x + d) : (x + ((d > 0.0) ? tol1 : -tol1));
    double const fu = objective(u);

    if (fu <= fx) {
      if (u < x) b = x; else a = x;
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) {
        v = w; fv = fw;
        w = u; fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u; fv = fu;
      }
    }
  }
  return x;
}

// Central-difference numerical gradient of a scalar objective, matching the
// numerical-differentiation fallback `stats::optim(method = "L-BFGS-B")` uses
// when called without an analytical `gr` (the former call site in
// `Fastcpd::GetOptimizedCostResult` never supplied one -- `cost_gradient_` /
// `cost_hessian_` were used only as the *trigger* for entering this
// gradient/Hessian-driven warm-start path, not as the gradient source: that
// user-supplied `cost_gradient_` is the SEN per-observation gradient
// (evaluated at a single row, for the sequential-gradient-descent update) and
// is not the gradient of the full-segment objective `cost_function_sen_`
// sums over -- substituting it here would pair an inconsistent (f, grad f)
// with the optimizer). `h` mirrors `optim`'s default `ndeps`-scale step.
template <typename Objective>
arma::colvec NumericalGradient(Objective const& objective,
                               arma::colvec const& theta,
                               double const h = 1e-6) {
  int const n = static_cast<int>(theta.n_elem);
  arma::colvec gradient(n);
  arma::colvec theta_perturbed = theta;
  for (int i = 0; i < n; i++) {
    double const theta_i = theta(i);
    theta_perturbed(i) = theta_i + h;
    double const f_plus = objective(theta_perturbed);
    theta_perturbed(i) = theta_i - h;
    double const f_minus = objective(theta_perturbed);
    theta_perturbed(i) = theta_i;
    gradient(i) = (f_plus - f_minus) / (2.0 * h);
  }
  return gradient;
}

// Projected backtracking line search satisfying the Armijo / sufficient-
// decrease condition (Nocedal & Wright, "Numerical Optimization" 2nd ed.,
// Algorithm 3.1; Bertsekas, "On the Goldstein-Levitin-Polyak gradient
// projection method", IEEE Trans. Automatic Control 21(2), 1976) — the
// step-size rule for the bound-constrained BFGS iteration below. Always
// leaves *theta_new/*f_new/*grad_new set to the last trial point, so the
// caller can proceed even if no step satisfied the condition within the
// attempt budget. Trial points are clamped to the box [lower, upper]
// (Bertsekas 1976) rather than capping the step length to stay strictly
// inside it -- see `BfgsMinimize` for why clamping is required for the
// active-set logic to engage cleanly at either bound.
template <typename Objective>
void ArmijoLineSearch(Objective const& objective, arma::colvec const& theta,
                      double const f, arma::colvec const& grad,
                      arma::colvec const& direction, arma::colvec const& lower,
                      arma::colvec const& upper, arma::colvec* theta_new,
                      double* f_new, arma::colvec* grad_new) {
  constexpr double kSufficientDecrease = 1e-4;
  constexpr double kBacktrack = 0.5;
  constexpr int kMaxAttempts = 50;
  double const directional_derivative = arma::dot(grad, direction);

  double step = 1.0;
  for (int attempt = 0; attempt < kMaxAttempts; attempt++) {
    arma::colvec const candidate =
        arma::min(arma::max(theta + step * direction, lower), upper);
    auto [f_candidate, grad_candidate] = objective(candidate);
    *theta_new = candidate;
    *f_new = f_candidate;
    *grad_new = grad_candidate;
    if (f_candidate <=
        f + kSufficientDecrease * step * directional_derivative) {
      return;
    }
    step *= kBacktrack;
  }
}

// Bound-constrained BFGS quasi-Newton minimisation over an axis-aligned box
// [lower, upper] (each bound may be +-infinity for unconstrained coordinates)
// via the gradient-projection / active-set method for simple bound
// constraints (Bertsekas, "Nonlinear Programming" 2nd ed., sec. 2.3; Nocedal
// & Wright, "Numerical Optimization" 2nd ed., sec. 16.7): maintain an
// inverse-Hessian approximation, restrict both the gradient and the
// quasi-Newton direction -B*grad to the "free" subspace -- coordinates that
// are either away from both bounds, or whose gradient sign indicates moving
// away from the active bound would decrease the objective -- choose a step by
// (boundary-respecting) Armijo line search, and refresh B via the secant/BFGS
// update. Used both by `fastcpd::garch::FitGarch` (replacing the vendored
// PORT `dsumsl` solver, with `lower = 0`, `upper = +Inf`) and by
// `Fastcpd::GetOptimizedCostResult` (replacing `stats::optim(method =
// "L-BFGS-B")` for the gradient/Hessian-driven PELT warm start).
//
// A coordinate at a bound (theta(i) ~= lower(i) with grad(i) >= 0, or
// theta(i) ~= upper(i) with grad(i) <= 0) already satisfies the first-order
// (KKT) optimality condition for a boundary solution -- moving further into
// the bound cannot decrease the objective -- so it is frozen ("active") for
// this iteration. A plain unconstrained line search cannot converge at such a
// boundary solution -- the feasible step length collapses to zero as the
// iterate approaches the boundary, and the iteration stalls with a nonzero
// gradient.
template <typename Objective>
arma::colvec BfgsMinimize(Objective const& objective,
                          arma::colvec const& theta0, arma::colvec const& lower,
                          arma::colvec const& upper, int const maxit = 200,
                          double const grad_tol = 1e-8) {
  // Coordinates within kBoundTol of a bound are treated as "at the bound" --
  // the iterate approaches, but in floating point essentially never exactly
  // reaches, the boundary.
  constexpr double kBoundTol = 1e-10;
  int const n = static_cast<int>(theta0.n_elem);
  arma::mat const identity = arma::eye<arma::mat>(n, n);
  arma::mat inv_hessian = identity;

  arma::colvec theta = arma::min(arma::max(theta0, lower), upper);
  auto [f, grad] = objective(theta);
  arma::uvec previous_active_idx;
  bool has_previous_active_idx = false;

  for (int iter = 0; iter < maxit; iter++) {
    arma::uvec is_active(n, arma::fill::zeros);
    for (int i = 0; i < n; i++) {
      bool const at_lower =
          theta(i) - lower(i) <= kBoundTol && grad(i) >= 0.0;
      bool const at_upper =
          upper(i) - theta(i) <= kBoundTol && grad(i) <= 0.0;
      if (at_lower || at_upper) is_active(i) = 1;
    }
    arma::uvec const active_idx = arma::find(is_active);

    // The curvature pairs (s, y) accumulated into inv_hessian reflect movement
    // restricted to the *previous* working subspace (the complement of the
    // active set); once the working set changes, that information no longer
    // describes the curvature of the newly-freed/frozen coordinates and can
    // produce a direction that drives a just-freed, near-boundary coordinate
    // straight back through the bound — the fraction-to-the-boundary rule then
    // shrinks the step geometrically every iteration without ever revising the
    // active set (the coordinate never quite reaches the bound). Restarting
    // the quasi-Newton approximation from the identity on every active-set
    // change (standard practice for active-set quasi-Newton methods, e.g.
    // Nocedal & Wright sec. 16.7) discards the stale curvature and lets the
    // steepest-descent-like first step re-orient correctly in the new
    // subspace.
    bool const active_set_changed =
        !has_previous_active_idx ||
        active_idx.n_elem != previous_active_idx.n_elem ||
        arma::any(active_idx != previous_active_idx);
    if (active_set_changed) inv_hessian = identity;
    previous_active_idx = active_idx;
    has_previous_active_idx = true;

    arma::colvec grad_free = grad;
    grad_free.elem(active_idx).zeros();
    if (arma::norm(grad_free, 2) < grad_tol) break;

    arma::colvec direction = -(inv_hessian * grad_free);
    direction.elem(active_idx).zeros();

    arma::colvec theta_new;
    double f_new;
    arma::colvec grad_new;
    ArmijoLineSearch(objective, theta, f, grad, direction, lower, upper,
                     &theta_new, &f_new, &grad_new);

    arma::colvec const s = theta_new - theta;
    arma::colvec const y_diff = grad_new - grad;
    double const sy = arma::dot(s, y_diff);
    if (sy > 1e-12) {
      double const rho = 1.0 / sy;
      inv_hessian = (identity - rho * s * y_diff.t()) * inv_hessian *
                        (identity - rho * y_diff * s.t()) +
                    rho * s * s.t();
    }

    theta = theta_new;
    f = f_new;
    grad = grad_new;

    // Near a flat optimum the objective can be constant to machine precision
    // along the descent direction (f_candidate == f for every representable
    // step length), so the strict Armijo sufficient-decrease condition
    // f_candidate <= f + c*step*directional_derivative — whose right-hand side
    // is always *strictly* less than f for a descent direction — can never be
    // satisfied; ArmijoLineSearch then exhausts its attempt budget and returns
    // a point bit-identical to theta, and the outer loop would otherwise spin
    // forever reproducing the same state. A step that moved theta by less than
    // machine epsilon cannot make further measurable progress, so stop here —
    // the standard "no progress" companion to the gradient-norm criterion.
    if (arma::norm(s, 2) < 1e-14) break;
  }

  return theta;
}

}  // namespace fastcpd::optim

#endif  // FASTCPD_OPTIM_H_
