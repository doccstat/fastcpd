#ifndef FASTCPD_FAMILIES_LASSO_H_
#define FASTCPD_FAMILIES_LASSO_H_

#include "fastcpd_family.h"

namespace fastcpd::families {

struct LassoFamily : BaseFamily {
  static constexpr const char* name = "lasso";
  static constexpr bool has_optimized_run = false;
  static constexpr bool has_error_std_dev = true;
  static constexpr bool is_lasso = true;

  static unsigned int GetDataNDims(arma::mat const& data) {
    return data.n_cols - 1;
  }

  template <typename Solver>
  static void CreateSenParameters(Solver* solver) {
    solver->coefficients_.col(0) = solver->segment_coefficients_.row(0).t();
    solver->coefficients_sum_.col(0) = solver->segment_coefficients_.row(0).t();
    solver->hessian_.slice(0) =
        TailOuter(solver, 0) +
        solver->epsilon_in_hessian_ *
            arma::eye<arma::mat>(solver->parameters_count_,
                                 solver->parameters_count_);
  }

  template <typename Solver>
  static void UpdateSenParameters(Solver* solver) {
    unsigned int const segment_index = solver->segment_index_;
    arma::colvec coef_add = solver->segment_coefficients_.row(segment_index).t();
    arma::mat hessian_new =
        TailOuter(solver, solver->t - 1) +
        solver->epsilon_in_hessian_ *
            arma::eye<arma::mat>(solver->parameters_count_,
                                 solver->parameters_count_);

    std::memcpy(solver->coefficients_.colptr(solver->pruned_set_size_ - 1),
                coef_add.memptr(), sizeof(double) * solver->parameters_count_);
    std::memcpy(solver->coefficients_sum_.colptr(solver->pruned_set_size_ - 1),
                coef_add.memptr(), sizeof(double) * solver->parameters_count_);
    std::memcpy(solver->hessian_.slice(solver->pruned_set_size_ - 1).memptr(),
                hessian_new.memptr(),
                sizeof(double) * solver->parameters_count_ *
                    solver->parameters_count_);
  }

  template <typename Solver>
  static void GetNllPelt(Solver* solver, unsigned int const segment_start,
                         unsigned int const segment_end, bool const cv,
                         std::optional<arma::colvec> const& start) {
    unsigned int const p = solver->data_.n_cols - 1;
    if (segment_start == segment_end) {
      solver->result_coefficients_ = arma::zeros<arma::colvec>(p);
      solver->result_residuals_ = cv ? arma::mat() : arma::zeros<arma::mat>(1, 1);
      solver->result_value_ = 0.0;
      return;
    }
    arma::mat const seg = solver->data_.rows(segment_start, segment_end);
    arma::colvec const y = seg.col(0);
    arma::mat const X = seg.cols(1, seg.n_cols - 1);
    if (cv) {
      auto const [beta, rss] = FitLassoCV(X, y);
      solver->result_coefficients_ = beta;
      solver->result_residuals_ = arma::mat();
      // Match original cv=true behavior: report RSS (not halved)
      solver->result_value_ = rss;
    } else {
      double const lambda =
          solver->lasso_penalty_base_ /
          std::sqrt(static_cast<double>(segment_end - segment_start + 1));
      auto const [beta, rss] = FitLassoFixed(X, y, lambda);
      solver->result_coefficients_ = beta;
      solver->result_residuals_ = y - X * beta;
      solver->result_value_ = rss / 2.0;
    }
  }

  template <typename Solver>
  static void GetNllPeltValue(Solver* solver, unsigned int const segment_start,
                              unsigned int const segment_end, bool const cv,
                              std::optional<arma::colvec> const& start) {
    GetNllPelt(solver, segment_start, segment_end, cv, start);
  }

  template <typename Solver>
  static double GetNllSen(Solver* solver, unsigned int const segment_start,
                          unsigned int const segment_end,
                          arma::colvec const& theta) {
    double cost = 0.0;
    for (unsigned int row = segment_start; row <= segment_end; row++) {
      double const residual = solver->data_(row, 0) - DotTail(solver, row, theta);
      cost += residual * residual;
    }
    return cost / 2.0 +
           solver->lasso_penalty_base_ /
               std::sqrt(segment_end - segment_start + 1) *
               arma::accu(arma::abs(theta));
  }

  template <typename Solver>
  static arma::colvec GetGradient(Solver* solver,
                                  unsigned int const segment_start,
                                  unsigned int const segment_end,
                                  arma::colvec const& theta) {
    double const y = solver->data_(segment_end, 0);
    double const eta = DotTail(solver, segment_end, theta);
    return ScaledTail(solver, segment_end, eta - y);
  }

  template <typename Solver>
  static arma::mat GetHessian(Solver* solver, unsigned int const segment_start,
                              unsigned int const segment_end,
                              arma::colvec const& theta) {
    return TailOuter(solver, segment_end);
  }

  // Soft-threshold: sign(z) * max(|z| - t, 0).
  static double SoftThreshold(double const z, double const t) {
    return std::copysign(std::max(std::abs(z) - t, 0.0), z);
  }

  // Coordinate descent lasso at fixed lambda.
  // Minimizes (1/(2n))||y - Xβ||² + lambda * ||β||₁.
  // Returns {beta, rss} where rss = ||y - Xβ||².
  static std::pair<arma::colvec, double> FitLassoFixed(
      arma::mat const& X, arma::colvec const& y, double const lambda) {
    unsigned int const n = X.n_rows, p = X.n_cols;
    arma::colvec beta = arma::zeros<arma::colvec>(p);
    arma::colvec residual = y;
    arma::rowvec const col_sq = arma::sum(arma::square(X), 0);
    for (unsigned int iter = 0; iter < 1000; iter++) {
      double max_change = 0.0;
      for (unsigned int j = 0; j < p; j++) {
        double const csq = col_sq(j);
        if (csq < 1e-14) continue;
        double const rho =
            (arma::dot(X.col(j), residual) + csq * beta(j)) / n;
        double const beta_new = SoftThreshold(rho, lambda) * n / csq;
        double const delta = beta_new - beta(j);
        if (std::abs(delta) > 1e-14) {
          residual -= X.col(j) * delta;
          beta(j) = beta_new;
          if (std::abs(delta) > max_change) max_change = std::abs(delta);
        }
      }
      if (max_change < 1e-7) break;
    }
    return {beta, arma::dot(residual, residual)};
  }

  // Lasso with 5-fold CV to select lambda (analogous to cv.glmnet lambda.1se).
  // Returns {beta, rss} at the selected lambda.
  static std::pair<arma::colvec, double> FitLassoCV(
      arma::mat const& X, arma::colvec const& y) {
    unsigned int const n = X.n_rows;
    double const lambda_max = arma::max(arma::abs(X.t() * y)) / n;
    if (lambda_max <= 0.0) return FitLassoFixed(X, y, 0.0);

    unsigned int constexpr kNLambda = 50;
    unsigned int const k_folds = std::min(5u, n);
    arma::colvec const lambdas = arma::exp(arma::linspace(
        std::log(lambda_max), std::log(lambda_max / 1000.0), kNLambda));

    // Deterministic fold assignment: sequential cyclic (balanced, no randomness).
    // Randomising the fold order with arma::shuffle introduces non-reproducibility
    // because Armadillo's RNG state depends on prior C++ calls in the same
    // session, making results vary with test execution order.
    arma::uvec fold_id(n);
    for (unsigned int i = 0; i < n; i++) fold_id(i) = i % k_folds;

    arma::mat fold_errors(k_folds, kNLambda, arma::fill::zeros);
    for (unsigned int fold = 0; fold < k_folds; fold++) {
      arma::uvec const tr = arma::find(fold_id != fold);
      arma::uvec const te = arma::find(fold_id == fold);
      if (tr.n_elem < 2 || te.is_empty()) continue;
      arma::mat const Xtr = X.rows(tr), Xte = X.rows(te);
      arma::colvec const ytr = y(tr), yte = y(te);
      for (unsigned int l = 0; l < kNLambda; l++) {
        arma::colvec const b = FitLassoFixed(Xtr, ytr, lambdas(l)).first;
        arma::colvec const resid = yte - Xte * b;
        fold_errors(fold, l) = arma::dot(resid, resid) / te.n_elem;
      }
    }

    arma::rowvec const mean_err = arma::mean(fold_errors, 0);
    arma::rowvec const se_err =
        arma::stddev(fold_errors, 0) / std::sqrt(static_cast<double>(k_folds));

    arma::uword const idx_min = mean_err.index_min();
    double const threshold = mean_err(idx_min) + se_err(idx_min);
    // lambda.1se: largest lambda with CV error within 1 SE of minimum.
    // Lambdas are decreasing, so scan from left (largest) to find first ≤ threshold.
    arma::uword idx_1se = idx_min;
    for (unsigned int l = 0; l < kNLambda; l++) {
      if (mean_err(l) <= threshold) { idx_1se = l; break; }
    }

    return FitLassoFixed(X, y, lambdas(idx_1se));
  }
};

}  // namespace fastcpd::families

#endif  // FASTCPD_FAMILIES_LASSO_H_
