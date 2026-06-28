#ifndef FASTCPD_FAMILIES_ARMA_H_
#define FASTCPD_FAMILIES_ARMA_H_

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "fastcpd_family.h"

namespace fastcpd {
namespace families {

class ArmaFamily : public BaseFamily {
 public:
  static constexpr bool has_optimized_run = false;

  template <typename T>
  static void GetNllPeltValue(T* instance, unsigned int const segment_start,
                              unsigned int const segment_end, bool const cv,
                              std::optional<arma::colvec> const& start) {
    GetNllPelt(instance, segment_start, segment_end, cv, start);
  }

  template <typename T>
  static void GetNllPelt(T* instance, unsigned int const segment_start,
                         unsigned int const segment_end, bool const cv,
                         std::optional<arma::colvec> const& start) {
    arma::colvec const series =
        instance->data_.rows(segment_start, segment_end).col(0);
    unsigned int const p = instance->order_(0);
    unsigned int const q = instance->order_(1);
    unsigned int const n_params = p + q + 1;
    unsigned int const seg_len = series.n_elem;

    if (q == 0) {
      // Pure AR: exact OLS conditional MLE.
      if (seg_len <= p) {
        instance->result_value_ = 0.0;
        instance->result_coefficients_ = arma::zeros<arma::colvec>(n_params);
        instance->result_residuals_ = arma::mat(series);
        return;
      }
      arma::mat X(seg_len - p, p);
      for (unsigned int j = 0; j < p; j++)
        X.col(j) = series.rows(p - j - 1, seg_len - j - 2);
      arma::colvec const y = series.rows(p, seg_len - 1);
      arma::colvec phi;
      if (!arma::solve(phi, X.t() * X, X.t() * y, arma::solve_opts::no_approx))
        phi = arma::zeros<arma::colvec>(p);
      arma::colvec const res = y - X * phi;
      double const n_eff = static_cast<double>(seg_len - p);
      double const sigma2 = std::max(arma::dot(res, res) / n_eff, 1e-10);
      instance->result_value_ =
          n_eff / 2.0 * (std::log(2.0 * M_PI) + std::log(sigma2) + 1.0);
      instance->result_coefficients_ = arma::zeros<arma::colvec>(n_params);
      instance->result_coefficients_.rows(0, p - 1) = phi;
      instance->result_coefficients_(n_params - 1) = sigma2;
      instance->result_residuals_ = arma::mat(res);
      return;
    }

    unsigned int const m = std::max(p, q);
    if (seg_len <= m) {
      instance->result_value_ = 0.0;
      instance->result_coefficients_ = arma::zeros<arma::colvec>(n_params);
      instance->result_residuals_ = arma::mat(series);
      return;
    }

    // Levinson-Durbin: partial autocorrelations kappa in (-1,1) -> AR/MA coefficients.
    // Ensures stationarity (AR) / invertibility (MA) by construction.
    auto pac_to_ar = [](arma::colvec const& kappa) -> arma::colvec {
      unsigned int const n = kappa.n_elem;
      if (n == 0) return arma::colvec();
      arma::colvec phi = arma::zeros<arma::colvec>(n);
      phi(0) = kappa(0);
      for (unsigned int k = 1; k < n; k++) {
        arma::colvec const prev = phi.rows(0, k - 1);
        for (unsigned int j = 0; j < k; j++)
          phi(j) = prev(j) - kappa(k) * prev(k - 1 - j);
        phi(k) = kappa(k);
      }
      return phi;
    };

    // Inverse Levinson-Durbin: AR/MA coefficients -> partial autocorrelations.
    auto ar_to_pac = [](arma::colvec phi) -> arma::colvec {
      unsigned int const n = phi.n_elem;
      arma::colvec kappa = arma::zeros<arma::colvec>(n);
      for (unsigned int k = n; k-- > 0; ) {
        kappa(k) = phi(k);
        if (k == 0) break;
        double const d = 1.0 - kappa(k) * kappa(k);
        if (d < 1e-10) break;
        arma::colvec prev = arma::zeros<arma::colvec>(k);
        for (unsigned int j = 0; j < k; j++)
          prev(j) = (phi(j) + kappa(k) * phi(k - 1 - j)) / d;
        phi.rows(0, k - 1) = prev;
      }
      return kappa;
    };

    unsigned int const r = std::max(p, q + 1u);

    // Build companion matrix and MA noise vector from AR/MA coefficients.
    // Correct state-space form for ARMA: first COLUMN = AR coefs, superdiagonal = 1.
    // This is the form used by R's stats::arima (Jones 1980).
    auto make_fm_rv = [&](arma::colvec const& phi,
                          arma::colvec const& psi)
        -> std::pair<arma::mat, arma::colvec> {
      arma::mat Fm = arma::zeros<arma::mat>(r, r);
      for (unsigned int i = 0; i < p; i++) Fm(i, 0) = phi(i);
      for (unsigned int i = 0; i < r - 1; i++) Fm(i, i + 1) = 1.0;
      arma::colvec Rv = arma::zeros<arma::colvec>(r);
      Rv(0) = 1.0;
      for (unsigned int j = 0; j < q; j++) Rv(j + 1) = psi(j);
      return {Fm, Rv};
    };

    // Stationary initialization: solve P0 = Fm * P0 * Fm' + Rv * Rv'
    // via vec(P0) = (Fm kron Fm) vec(P0) + vec(Rv*Rv').
    auto stationary_cov = [&](arma::mat const& Fm,
                               arma::colvec const& Rv) -> arma::mat {
      unsigned int const rr = r * r;
      arma::mat const I_FkF =
          arma::eye<arma::mat>(rr, rr) - arma::kron(Fm, Fm);
      arma::colvec vec_P0;
      if (!arma::solve(vec_P0, I_FkF, arma::vectorise(Rv * Rv.t()),
                       arma::solve_opts::no_approx))
        return arma::eye<arma::mat>(r, r);
      return arma::reshape(vec_P0, r, r);
    };

    // Diffuse Kalman ML in unconstrained eta space.
    // AR: eta -> kappa=tanh(eta) -> phi=pac_to_ar(kappa) (stationarity)
    // MA: eta -> kappa=tanh(eta) -> psi=-pac_to_ar(kappa) (invertibility)
    auto kalman_nll = [&](arma::colvec const& eta) -> double {
      arma::colvec kap_phi(p), kap_psi(q);
      for (unsigned int i = 0; i < p; i++) kap_phi(i) = std::tanh(eta(i));
      for (unsigned int j = 0; j < q; j++) kap_psi(j) = std::tanh(eta(p + j));
      arma::colvec const phi = p > 0 ? pac_to_ar(kap_phi) : arma::colvec();
      arma::colvec const psi = -pac_to_ar(kap_psi);

      auto [Fm, Rv] = make_fm_rv(phi, psi);
      arma::mat P = stationary_cov(Fm, Rv);
      if (!P.is_finite()) return 1e10;

      arma::colvec a = arma::zeros<arma::colvec>(r);
      double sv2f = 0.0, slF = 0.0;
      for (unsigned int t = 0; t < seg_len; t++) {
        double const v = series(t) - a(0);
        double const F = P(0, 0);
        if (F <= 0.0 || !std::isfinite(F)) return 1e10;
        sv2f += v * v / F;
        slF += std::log(F);
        arma::colvec const Kv = Fm * (P.col(0) / F);
        a = Fm * a + Kv * v;
        P = Fm * (P - P.col(0) * P.row(0) / F) * Fm.t() + Rv * Rv.t();
      }
      double const s2 = sv2f / seg_len;
      if (s2 <= 0.0 || !std::isfinite(s2)) return 1e10;
      return 0.5 * seg_len * (1.0 + std::log(2.0 * M_PI) + std::log(s2))
             + 0.5 * slF;
    };

    // Hannan-Rissanen starting point in natural phi/psi space.
    unsigned int const p_approx = std::min(
        std::max(p + q + 5u, 2u * (p + q)),
        (seg_len > 2u ? (seg_len - 1u) / 2u : 1u));
    arma::colvec eps_init = arma::zeros<arma::colvec>(seg_len);
    if (p_approx >= 1 && seg_len > p_approx + 1) {
      arma::mat X_ar(seg_len - p_approx, p_approx);
      for (unsigned int j = 0; j < p_approx; j++)
        X_ar.col(j) = series.rows(p_approx - j - 1, seg_len - j - 2);
      arma::colvec phi_ar;
      if (arma::solve(phi_ar, X_ar.t() * X_ar,
                      X_ar.t() * series.rows(p_approx, seg_len - 1),
                      arma::solve_opts::no_approx)) {
        for (unsigned int t = p_approx; t < seg_len; t++) {
          eps_init(t) = series(t);
          for (unsigned int j = 0; j < p_approx; j++)
            eps_init(t) -= phi_ar(j) * series(t - j - 1);
        }
      }
    }
    arma::colvec params_nat = arma::zeros<arma::colvec>(p + q);
    if (seg_len > m + p + q) {
      arma::mat X_arma(seg_len - m, p + q);
      for (unsigned int j = 0; j < p; j++)
        X_arma.col(j) = series.rows(m - j - 1, seg_len - j - 2);
      for (unsigned int k = 0; k < q; k++)
        X_arma.col(p + k) = eps_init.rows(m - k - 1, seg_len - k - 2);
      arma::colvec params_hr;
      if (arma::solve(params_hr, X_arma.t() * X_arma,
                      X_arma.t() * series.rows(m, seg_len - 1),
                      arma::solve_opts::no_approx))
        params_nat = arma::clamp(params_hr, -0.97, 0.97);
    }

    // Convert HR natural-space starting point to unconstrained eta via PAC.
    auto safe_atanh = [](double x) -> double {
      return std::atanh(std::max(-0.97, std::min(0.97, x)));
    };
    arma::colvec eta0 = arma::zeros<arma::colvec>(p + q);
    if (p > 0) {
      arma::colvec const kap = ar_to_pac(params_nat.rows(0, p - 1));
      for (unsigned int i = 0; i < p; i++) eta0(i) = safe_atanh(kap(i));
    }
    {
      arma::colvec const kap = ar_to_pac(-params_nat.rows(p, p + q - 1));
      for (unsigned int j = 0; j < q; j++) eta0(p + j) = safe_atanh(kap(j));
    }

    // Nelder-Mead in unconstrained eta space (alpha=1, gamma=2, rho=0.5, sigma=0.5).
    {
      unsigned int const d = p + q;
      unsigned int const nv = d + 1;
      std::vector<arma::colvec> V(nv);
      std::vector<double> fV(nv);
      V[0] = eta0; fV[0] = kalman_nll(eta0);
      for (unsigned int i = 1; i <= d; i++) {
        V[i] = eta0;
        V[i](i - 1) += (std::abs(eta0(i - 1)) < 1e-4 ? 0.05
                                                        : 0.05 * std::abs(eta0(i - 1)));
        fV[i] = kalman_nll(V[i]);
      }
      std::vector<unsigned int> ord(nv);
      for (unsigned int it = 0; it < 500; it++) {
        std::iota(ord.begin(), ord.end(), 0);
        std::sort(ord.begin(), ord.end(),
                  [&](unsigned int a, unsigned int b) { return fV[a] < fV[b]; });
        if (fV[ord[nv - 1]] - fV[ord[0]] < 1e-8) break;
        arma::colvec c = arma::zeros<arma::colvec>(d);
        for (unsigned int i = 0; i < d; i++) c += V[ord[i]];
        c /= static_cast<double>(d);
        arma::colvec xr = 2.0 * c - V[ord[nv - 1]]; double fr = kalman_nll(xr);
        if (fr < fV[ord[0]]) {
          arma::colvec xe = 3.0 * c - 2.0 * V[ord[nv - 1]]; double fe = kalman_nll(xe);
          if (fe < fr) { V[ord[nv-1]] = xe; fV[ord[nv-1]] = fe; }
          else          { V[ord[nv-1]] = xr; fV[ord[nv-1]] = fr; }
        } else if (fr < fV[ord[d - 1]]) {
          V[ord[nv - 1]] = xr; fV[ord[nv - 1]] = fr;
        } else {
          bool shrink = false;
          if (fr < fV[ord[nv - 1]]) {
            arma::colvec xc = 0.5 * (c + xr); double fc = kalman_nll(xc);
            if (fc <= fr) { V[ord[nv-1]] = xc; fV[ord[nv-1]] = fc; } else shrink = true;
          } else {
            arma::colvec xc = 0.5 * (c + V[ord[nv-1]]); double fc = kalman_nll(xc);
            if (fc < fV[ord[nv-1]]) { V[ord[nv-1]] = xc; fV[ord[nv-1]] = fc; } else shrink = true;
          }
          if (shrink) {
            for (unsigned int i = 1; i < nv; i++) {
              V[ord[i]] = 0.5 * (V[ord[0]] + V[ord[i]]);
              fV[ord[i]] = kalman_nll(V[ord[i]]);
            }
          }
        }
      }
      unsigned int best = 0;
      for (unsigned int i = 1; i < nv; i++) if (fV[i] < fV[best]) best = i;
      eta0 = V[best];
    }

    // Decode optimal eta -> phi*, psi*.
    arma::colvec phi_star;
    arma::colvec psi_star(q);
    {
      if (p > 0) {
        arma::colvec kap_phi(p);
        for (unsigned int i = 0; i < p; i++) kap_phi(i) = std::tanh(eta0(i));
        phi_star = pac_to_ar(kap_phi);
      }
      arma::colvec kap_psi(q);
      for (unsigned int j = 0; j < q; j++) kap_psi(j) = std::tanh(eta0(p + j));
      psi_star = -pac_to_ar(kap_psi);
    }

    // Final KF pass with stationary init: recover sigma2 and Kalman innovations.
    {
      auto [Fm, Rv] = make_fm_rv(phi_star, psi_star);
      arma::mat P = stationary_cov(Fm, Rv);
      arma::colvec a = arma::zeros<arma::colvec>(r);
      arma::colvec innov = arma::zeros<arma::colvec>(seg_len);
      double sv2f = 0.0;
      for (unsigned int t = 0; t < seg_len; t++) {
        double const v = series(t) - a(0);
        double const F = P(0, 0);
        innov(t) = v;
        if (F <= 0.0 || !std::isfinite(F)) break;
        sv2f += v * v / F;
        arma::colvec const Kv = Fm * (P.col(0) / F);
        a = Fm * a + Kv * v;
        P = Fm * (P - P.col(0) * P.row(0) / F) * Fm.t() + Rv * Rv.t();
      }
      double const sigma2 = std::max(sv2f / seg_len, 1e-10);
      arma::colvec theta_full = arma::zeros<arma::colvec>(n_params);
      if (p > 0) theta_full.rows(0, p - 1) = phi_star;
      theta_full.rows(p, p + q - 1) = psi_star;
      theta_full(n_params - 1) = sigma2;
      instance->result_value_ = kalman_nll(eta0);
      instance->result_coefficients_ = theta_full;
      instance->result_residuals_ = arma::mat(innov);
    }
  }

  template <typename T>
  static double GetNllSen(T* instance, unsigned int const segment_start,
                          unsigned int const segment_end,
                          arma::colvec const& theta) {
    arma::mat data_segment = instance->data_.rows(segment_start, segment_end);
    arma::colvec reversed_theta = arma::reverse(theta);
    if (data_segment.n_rows < arma::max(instance->order_) + 1) {
      return 0;
    }
    arma::colvec variance_term = arma::zeros(data_segment.n_rows);
    for (unsigned int i = arma::max(instance->order_); i < data_segment.n_rows; i++) {
      variance_term(i) = data_segment(i, 0) -
                         arma::dot(reversed_theta.rows(instance->order_(1) + 1, arma::sum(instance->order_)),
                                   data_segment.rows(i - instance->order_(0), i - 1)) -
                         arma::dot(reversed_theta.rows(1, instance->order_(1)),
                                   variance_term.rows(i - instance->order_(1), i - 1));
    }
    return (std::log(2.0 * M_PI) + std::log(theta(arma::sum(instance->order_)))) *
               (data_segment.n_rows - 2) / 2.0 +
           arma::dot(variance_term, variance_term) / 2.0 / theta(arma::sum(instance->order_));
  }

  template <typename T>
  static arma::colvec GetGradient(T* instance, unsigned int const segment_start,
                                  unsigned int const segment_end,
                                  arma::colvec const& theta) {
    arma::mat const data_segment = instance->data_.rows(segment_start, segment_end);
    unsigned int const segment_length = segment_end - segment_start + 1;
    arma::mat reversed_data = arma::reverse(data_segment, 0);
    arma::colvec reversed_theta = arma::reverse(theta);
    if (segment_length < arma::max(instance->order_) + 1) {
      return arma::ones(theta.n_elem);
    }
    arma::colvec variance_term = arma::zeros(segment_length);
    for (unsigned int i = arma::max(instance->order_); i < segment_length; i++) {
      variance_term(i) = data_segment(i, 0) -
                         arma::dot(reversed_theta.rows(instance->order_(1) + 1, arma::sum(instance->order_)),
                                   data_segment.rows(i - instance->order_(0), i - 1)) -
                         arma::dot(reversed_theta.rows(1, instance->order_(1)),
                                   variance_term.rows(i - instance->order_(1), i - 1));
    }
    arma::colvec reversed_variance_term = arma::reverse(variance_term);
    arma::mat phi_coefficient = arma::zeros(segment_length, instance->order_(0)),
              psi_coefficient = arma::zeros(segment_length, instance->order_(1));
    for (unsigned int i = arma::max(instance->order_); i < segment_length; i++) {
      phi_coefficient.row(i) =
          -reversed_data
               .rows(segment_length - i, segment_length - i + instance->order_(0) - 1)
               .t() -
          reversed_theta.rows(1, instance->order_(1)).t() *
              phi_coefficient.rows(i - instance->order_(1), i - 1);
    }
    for (unsigned int i = instance->order_(1); i < segment_length; i++) {
      psi_coefficient.row(i) =
          -reversed_variance_term
               .rows(segment_length - i, segment_length - i + instance->order_(1) - 1)
               .t() -
          reversed_theta.rows(1, instance->order_(1)).t() *
              psi_coefficient.rows(i - instance->order_(1), i - 1);
    }
    arma::colvec gradient = arma::zeros(arma::sum(instance->order_) + 1);
    gradient.rows(0, instance->order_(0) - 1) =
        phi_coefficient.row(segment_length - 1).t() *
        variance_term(segment_length - 1) / theta(arma::sum(instance->order_));
    gradient.rows(instance->order_(0), arma::sum(instance->order_) - 1) =
        psi_coefficient.row(segment_length - 1).t() *
        variance_term(segment_length - 1) / theta(arma::sum(instance->order_));
    gradient(arma::sum(instance->order_)) = 1.0 / 2.0 / theta(arma::sum(instance->order_)) -
                            pow(variance_term(segment_length - 1), 2) / 2.0 /
                                pow(theta(arma::sum(instance->order_)), 2);
    return gradient;
  }

  template <typename T>
  static arma::mat GetHessian(T* instance, unsigned int const segment_start,
                              unsigned int const segment_end,
                              arma::colvec const& theta) {
    arma::mat const data_segment = instance->data_.rows(segment_start, segment_end);
    unsigned int const segment_length = segment_end - segment_start + 1;
    arma::mat reversed_data = arma::reverse(data_segment, 0);
    arma::colvec reversed_theta = arma::reverse(theta);
    if (segment_length < arma::max(instance->order_) + 1) {
      return arma::eye(theta.n_elem, theta.n_elem);
    }
    arma::colvec variance_term = arma::zeros(segment_length);
    for (unsigned int i = arma::max(instance->order_); i < segment_length; i++) {
      variance_term(i) = data_segment(i, 0) -
                         arma::dot(reversed_theta.rows(instance->order_(1) + 1, arma::sum(instance->order_)),
                                   data_segment.rows(i - instance->order_(0), i - 1)) -
                         arma::dot(reversed_theta.rows(1, instance->order_(1)),
                                   variance_term.rows(i - instance->order_(1), i - 1));
    }
    arma::colvec reversed_variance_term = arma::reverse(variance_term);
    arma::mat phi_coefficient = arma::zeros(segment_length, instance->order_(0)),
              psi_coefficient = arma::zeros(segment_length, instance->order_(1));
    for (unsigned int i = arma::max(instance->order_); i < segment_length; i++) {
      phi_coefficient.row(i) =
          -reversed_data
               .rows(segment_length - i, segment_length - i + instance->order_(0) - 1)
               .t() -
          reversed_theta.rows(1, instance->order_(1)).t() *
              phi_coefficient.rows(i - instance->order_(1), i - 1);
    }
    for (unsigned int i = instance->order_(1); i < segment_length; i++) {
      psi_coefficient.row(i) =
          -reversed_variance_term
               .rows(segment_length - i, segment_length - i + instance->order_(1) - 1)
               .t() -
          reversed_theta.rows(1, instance->order_(1)).t() *
              psi_coefficient.rows(i - instance->order_(1), i - 1);
    }
    arma::mat reversed_coef_phi = arma::reverse(phi_coefficient, 0),
              reversed_coef_psi = arma::reverse(psi_coefficient, 0);
    arma::cube phi_psi_coefficient = arma::zeros(instance->order_(1), instance->order_(0), segment_length),
               psi_psi_coefficient = arma::zeros(instance->order_(1), instance->order_(1), segment_length);
    for (unsigned int i = instance->order_(1); i < segment_length; i++) {
      arma::mat phi_psi_coefficient_part = arma::zeros(instance->order_(1), instance->order_(0)),
                psi_psi_coefficient_part = arma::zeros(instance->order_(1), instance->order_(1));
      for (unsigned int j = 1; j <= instance->order_(1); j++) {
        phi_psi_coefficient_part +=
            phi_psi_coefficient.slice(i - j) * theta(instance->order_(0) - 1 + j);
      }
      phi_psi_coefficient.slice(i) =
          -reversed_coef_phi.rows(segment_length - i,
                                  segment_length - i + instance->order_(1) - 1) -
          phi_psi_coefficient_part;
      for (unsigned int j = 1; j <= instance->order_(1); j++) {
        psi_psi_coefficient_part +=
            psi_psi_coefficient.slice(i - j) * theta(instance->order_(0) - 1 + j);
      }
      psi_psi_coefficient.slice(i) =
          -reversed_coef_psi.rows(segment_length - i,
                                  segment_length - i + instance->order_(1) - 1) -
          reversed_coef_psi
              .rows(segment_length - i, segment_length - i + instance->order_(1) - 1)
              .t() -
          psi_psi_coefficient_part;
    }
    arma::mat hessian = arma::zeros(arma::sum(instance->order_) + 1, arma::sum(instance->order_) + 1);
    hessian.submat(0, 0, instance->order_(0) - 1, instance->order_(0) - 1) =
        phi_coefficient.row(segment_length - 1).t() *
        phi_coefficient.row(segment_length - 1) / theta(arma::sum(instance->order_));
    hessian.submat(0, instance->order_(0), instance->order_(0) - 1, arma::sum(instance->order_) - 1) =
        (phi_psi_coefficient.slice(segment_length - 1).t() *
             variance_term(segment_length - 1) +
         phi_coefficient.row(segment_length - 1).t() *
             psi_coefficient.row(segment_length - 1)) /
        theta(arma::sum(instance->order_));
    hessian.submat(instance->order_(0), 0, arma::sum(instance->order_) - 1, instance->order_(0) - 1) =
        hessian.submat(0, instance->order_(0), instance->order_(0) - 1, arma::sum(instance->order_) - 1).t();
    hessian.submat(0, arma::sum(instance->order_), instance->order_(0) - 1, arma::sum(instance->order_)) =
        -phi_coefficient.row(segment_length - 1).t() *
        variance_term(segment_length - 1) / theta(arma::sum(instance->order_)) /
        theta(arma::sum(instance->order_));
    hessian.submat(arma::sum(instance->order_), 0, arma::sum(instance->order_), instance->order_(0) - 1) =
        hessian.submat(0, arma::sum(instance->order_), instance->order_(0) - 1, arma::sum(instance->order_)).t();
    hessian.submat(instance->order_(0), instance->order_(0), arma::sum(instance->order_) - 1, arma::sum(instance->order_) - 1) =
        (psi_coefficient.row(segment_length - 1).t() *
             psi_coefficient.row(segment_length - 1) +
         psi_psi_coefficient.slice(segment_length - 1) *
             variance_term(segment_length - 1)) /
        theta(arma::sum(instance->order_));
    hessian.submat(instance->order_(0), arma::sum(instance->order_), arma::sum(instance->order_) - 1, arma::sum(instance->order_)) =
        -psi_coefficient.row(segment_length - 1).t() *
        variance_term(segment_length - 1) / theta(arma::sum(instance->order_)) /
        theta(arma::sum(instance->order_));
    hessian.submat(arma::sum(instance->order_), instance->order_(0), arma::sum(instance->order_), arma::sum(instance->order_) - 1) =
        hessian.submat(instance->order_(0), arma::sum(instance->order_), arma::sum(instance->order_) - 1, arma::sum(instance->order_)).t();
    hessian(arma::sum(instance->order_), arma::sum(instance->order_)) =
        pow(variance_term(segment_length - 1), 2) / pow(theta(arma::sum(instance->order_)), 3) -
        1.0 / 2.0 / pow(theta(arma::sum(instance->order_)), 2);
    return hessian;
  }

  template <typename Solver>
  static void CreateSenParameters(Solver* solver) {
    solver->coefficients_.col(0) = solver->segment_coefficients_.row(0).t();
    solver->coefficients_sum_.col(0) = solver->segment_coefficients_.row(0).t();
    solver->hessian_.slice(0) =
        arma::zeros<arma::mat>(solver->parameters_count_, solver->parameters_count_);
  }

  template <typename Solver>
  static void UpdateSenParameters(Solver* solver) {
    int const segment_index =
        arma::index_max(arma::find(solver->segment_indices_ <= solver->t - 1));
    arma::colvec coef_add = solver->segment_coefficients_.row(segment_index).t();
    arma::mat hessian_new =
        ArmaFamily::GetHessian(solver, 0, solver->t - 1, coef_add);
    std::memcpy(solver->coefficients_.colptr(solver->pruned_set_size_ - 1),
                coef_add.memptr(),
                sizeof(double) * solver->parameters_count_);
    std::memcpy(solver->coefficients_sum_.colptr(solver->pruned_set_size_ - 1),
                coef_add.memptr(),
                sizeof(double) * solver->parameters_count_);
    std::memcpy(
        solver->hessian_.slice(solver->pruned_set_size_ - 1).memptr(),
        hessian_new.memptr(),
        sizeof(double) * solver->parameters_count_ * solver->parameters_count_);
  }
};

} // namespace families
} // namespace fastcpd

#endif
