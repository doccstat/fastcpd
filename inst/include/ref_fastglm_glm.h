// This is a modified copy of fastglm/inst/include/glm.h

#ifndef GLM_H
#define GLM_H

#include "ref_fastglm_glm_base.h"

using Eigen::ArrayXd;
using Eigen::BDCSVD;
using Eigen::ColPivHouseholderQR;
using Eigen::ComputeThinU;
using Eigen::ComputeThinV;
using Eigen::EigenBase;
using Eigen::FullPivHouseholderQR;
using Eigen::HouseholderQR;
using Eigen::JacobiSVD;
using Eigen::LDLT;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SelfAdjointEigenSolver;
using Eigen::SelfAdjointView;
using Eigen::TriangularView;
using Eigen::Upper;
using Eigen::VectorXd;

class glm : public GlmBase<Eigen::VectorXd,
                           Eigen::MatrixXd>  // Eigen::SparseVector<double>
{
 protected:
  typedef double Double;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Map<Matrix const> MapMat;
  typedef Eigen::Map<Vector const> MapVec;
  typedef Eigen::Ref<Matrix const> const ConstGenericMatrix;
  typedef Eigen::Ref<Vector const> const ConstGenericVector;
  typedef Eigen::SparseMatrix<double> SpMat;
  typedef Eigen::SparseVector<double> SparseVector;

  typedef MatrixXd::Index Index;
  typedef MatrixXd::Scalar Scalar;
  typedef MatrixXd::RealScalar RealScalar;
  typedef ColPivHouseholderQR<MatrixXd>::PermutationType Permutation;
  typedef Permutation::IndicesType Indices;

  Map<MatrixXd> const X;
  Map<VectorXd> const Y;
  Map<VectorXd> const weights;
  Map<VectorXd> const offset;

  NumericVector (*variance_fun)(Eigen::VectorXd const&);
  NumericVector (*mu_eta_fun)(Eigen::VectorXd const&);
  NumericVector (*linkinv)(Eigen::VectorXd const&);
  NumericVector (*dev_resids_fun)(Map<VectorXd> const&, Eigen::VectorXd const&,
                                  Map<VectorXd> const&);
  bool (*valideta)(Eigen::VectorXd const&);
  bool (*validmu)(Eigen::VectorXd const&);

  double tol;
  int maxit;
  int type;
  bool is_big_matrix;
  int rank;

  FullPivHouseholderQR<MatrixXd> FPQR;
  ColPivHouseholderQR<MatrixXd> PQR;
  BDCSVD<MatrixXd> bSVD;
  HouseholderQR<MatrixXd> QR;
  LLT<MatrixXd> Ch;
  LDLT<MatrixXd> ChD;
  JacobiSVD<MatrixXd> UDV;

  SelfAdjointEigenSolver<MatrixXd> eig;

  Permutation Pmat;
  MatrixXd Rinv;
  VectorXd effects;

  RealScalar threshold() const {
    // return m_usePrescribedThreshold ? m_prescribedThreshold
    //: numeric_limits<double>::epsilon() * nvars;
    return numeric_limits<double>::epsilon() * nvars;
  }

  // from RcppEigen
  inline ArrayXd Dplus(ArrayXd const& d) {
    ArrayXd di(d.size());
    double comp(d.maxCoeff() * threshold());
    for (int j = 0; j < d.size(); ++j) di[j] = (d[j] < comp) ? 0. : 1. / d[j];
    rank = (di != 0.).count();
    return di;
  }

  MatrixXd XtWX() const {
    if (!is_big_matrix || true) {
      return MatrixXd(nvars, nvars)
          .setZero()
          .selfadjointView<Lower>()
          .rankUpdate((w.asDiagonal() * X).adjoint());
    } else {
      MatrixXd retmat(MatrixXd::Zero(nvars, nvars));
      VectorXd wsquare = w.array().square();
      for (int j = 0; j < nvars; ++j) {
        for (int k = j; k < nvars; ++k) {
          for (int i = 0; i < nobs; ++i) {
            retmat(j, k) += X(i, j) * X(i, k) * wsquare(i);
          }
          retmat(k, j) = retmat(j, k);
        }
      }
      return retmat;
    }
  }

  virtual void update_mu_eta() {
    NumericVector mu_eta_nv = mu_eta_fun(eta);

    std::copy(mu_eta_nv.begin(), mu_eta_nv.end(), mu_eta.data());
  }

  virtual void update_var_mu() {
    NumericVector var_mu_nv = variance_fun(mu);

    std::copy(var_mu_nv.begin(), var_mu_nv.end(), var_mu.data());
  }

  virtual void update_mu() {
    // mu <- linkinv(eta <- eta + offset)
    NumericVector mu_nv = linkinv(eta);

    std::copy(mu_nv.begin(), mu_nv.end(), mu.data());
  }

  virtual void update_eta() {
    // eta <- drop(x %*% start)

    if (type == 0) {
      // VectorXd effects(PQR.householderQ().adjoint() * y);
      if (rank == nvars) {
        eta = X * beta + offset;
      } else {
        // eta = PQR.householderQ() * effects + offset;
        eta = X * beta + offset;
      }
    } else if (type == 1) {
      eta = X * beta + offset;
    } else if (type == 2) {
      eta = X * beta + offset;
    } else if (type == 3) {
      eta = X * beta + offset;
    } else if (type == 4) {
      if (rank == nvars) {
        eta = X * beta + offset;
      } else {
        // std::cout << FPQR.matrixQ().cols() << " " << effects.size() <<
        // std::endl; eta = FPQR.matrixQ() * effects + offset;
        eta = X * beta + offset;
      }
    } else {
      eta = X * beta + offset;
    }
  }

  virtual void update_z() {
    // z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
    z = (eta - offset).array() + (Y - mu).array() / mu_eta.array();
  }

  virtual void update_w() {
    // w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
    w = (weights.array() * mu_eta.array().square() / var_mu.array())
            .array()
            .sqrt();
  }

  virtual void update_dev_resids() {
    devold = dev;
    NumericVector dev_resids = dev_resids_fun(Y, mu, weights);
    dev = sum(dev_resids);
  }

  virtual void update_dev_resids_dont_update_old() {
    NumericVector dev_resids = dev_resids_fun(Y, mu, weights);
    dev = sum(dev_resids);
  }

  virtual void step_halve() {
    // take half step
    beta = 0.5 * (beta.array() + beta_prev.array());

    update_eta();

    update_mu();
  }

  virtual void run_step_halving(int& iterr) {
    // check for infinite deviance
    if (std::isinf(dev)) {
      int itrr = 0;
      while (std::isinf(dev)) {
        ++itrr;
        if (itrr > maxit) {
          break;
        }

        // std::cout << "half step (infinite)!" << itrr << std::endl;

        step_halve();

        // update deviance
        update_dev_resids_dont_update_old();
      }
    }

    // check for boundary violations
    if (!(valideta(eta) && validmu(mu))) {
      int itrr = 0;
      while (!(valideta(eta) && validmu(mu))) {
        ++itrr;
        if (itrr > maxit) {
          break;
        }

        // std::cout << "half step (boundary)!" << itrr << std::endl;

        step_halve();
      }

      update_dev_resids_dont_update_old();
    }

    // check for increasing deviance
    // std::abs(deviance - deviance_prev) / (0.1 + std::abs(deviance)) <
    // tol_irls
    if ((dev - devold) / (0.1 + std::abs(dev)) >= tol && iterr > 0) {
      int itrr = 0;

      // std::cout << "dev:" << deviance << "dev prev:" << deviance_prev <<
      // std::endl;

      while ((dev - devold) / (0.1 + std::abs(dev)) >= -tol) {
        ++itrr;
        if (itrr > maxit) {
          break;
        }

        // std::cout << "half step (increasing dev)!" << itrr << std::endl;

        step_halve();

        update_dev_resids_dont_update_old();
      }
    }
  }

  // much of solve_wls() comes directly
  // from the source code of the RcppEigen package
  virtual void solve_wls(int iter) {
    // lm ans(do_lm(X, Y, w, type));
    // wls ans(ColPivQR(X, z, w));

    // enum {ColPivQR_t = 0, QR_t, LLT_t, LDLT_t, SVD_t, SymmEigen_t, GESDD_t};

    beta_prev = beta;

    if (type == 0) {
      PQR.compute(w.asDiagonal() * X);  // decompose the model matrix
      Pmat = (PQR.colsPermutation());
      rank = PQR.rank();
      if (rank == nvars) {  // full rank case
        beta = PQR.solve((z.array() * w.array()).matrix());
        // m_fitted   = X * m_coef;
        // m_se       = Pmat * PQR.matrixQR().topRows(m_p).
        // triangularView<Upper>().solve(MatrixXd::Identity(nvars,
        // nvars)).rowwise().norm();
      } else {
        Rinv = (PQR.matrixQR()
                    .topLeftCorner(rank, rank)
                    .triangularView<Upper>()
                    .solve(MatrixXd::Identity(rank, rank)));
        effects =
            PQR.householderQ().adjoint() * (z.array() * w.array()).matrix();
        beta.head(rank) = Rinv * effects.head(rank);
        beta = Pmat * beta;

        // create fitted values from effects
        // (can't use X*m_coef if X is rank-deficient)
        effects.tail(nobs - rank).setZero();
        // m_fitted                          = PQR.householderQ() * effects;
        // m_se.head(m_r)                    = Rinv.rowwise().norm();
        // m_se                              = Pmat * m_se;
      }
    } else if (type == 1) {
      QR.compute(w.asDiagonal() * X);
      beta = QR.solve((z.array() * w.array()).matrix());
      // m_fitted                   = X * m_coef;
      // m_se                       = QR.matrixQR().topRows(m_p).
      // triangularView<Upper>().solve(I_p()).rowwise().norm();
    } else if (type == 2) {
      Ch.compute(XtWX().selfadjointView<Lower>());
      beta = Ch.solve((w.asDiagonal() * X).adjoint() *
                      (z.array() * w.array()).matrix());
      // m_fitted          = X * m_coef;
      // m_se              = Ch.matrixL().solve(I_p()).colwise().norm();
    } else if (type == 3) {
      ChD.compute(XtWX().selfadjointView<Lower>());
      Dplus(ChD.vectorD());  // to set the rank
      // FIXME: Check on the permutation in the LDLT and incorporate it in
      // the coefficients and the standard error computation.
      //	m_coef            = Ch.matrixL().adjoint().
      //	    solve(Dplus(D) * Ch.matrixL().solve(X.adjoint() * y));
      beta = ChD.solve((w.asDiagonal() * X).adjoint() *
                       (z.array() * w.array()).matrix());
      // m_fitted          = X * m_coef;
      // m_se              = Ch.solve(I_p()).diagonal().array().sqrt();
    } else if (type == 4) {
      FPQR.compute(w.asDiagonal() * X);  // decompose the model matrix
      Pmat = (FPQR.colsPermutation());
      rank = FPQR.rank();
      if (rank == nvars) {  // full rank case
        beta = FPQR.solve((z.array() * w.array()).matrix());
        // m_fitted   = X * m_coef;
        // m_se       = Pmat * PQR.matrixQR().topRows(m_p).
        // triangularView<Upper>().solve(MatrixXd::Identity(nvars,
        // nvars)).rowwise().norm();
      } else {
        Rinv = (FPQR.matrixQR()
                    .topLeftCorner(rank, rank)
                    .triangularView<Upper>()
                    .solve(MatrixXd::Identity(rank, rank)));
        effects = FPQR.matrixQ().adjoint() * (z.array() * w.array()).matrix();
        // std::cout << effects.transpose() << std::endl;
        beta.head(rank) = Rinv * effects.head(rank);
        beta = Pmat * beta;

        // create fitted values from effects
        // (can't use X*m_coef if X is rank-deficient)
        effects.tail(nobs - rank).setZero();
        // m_fitted                          = PQR.householderQ() * effects;
        // m_se.head(m_r)                    = Rinv.rowwise().norm();
        // m_se                              = Pmat * m_se;
      }
    } else if (type == 5) {
      bSVD.compute(w.asDiagonal() * X, ComputeThinU | ComputeThinV);

      rank = bSVD.rank();

      // if (rank == nvars)
      // {	// full rank case
      //     beta                     = bSVD.solve((z.array() *
      //     w.array()).matrix());
      // } else
      // {
      //
      // }

      beta = bSVD.solve((z.array() * w.array()).matrix());

      // FIXME: Check on the permutation in the LDLT and incorporate it in
      // the coefficients and the standard error computation.
      //	m_coef            = Ch.matrixL().adjoint().
      //	    solve(Dplus(D) * Ch.matrixL().solve(X.adjoint() * y));
      // m_fitted          = X * m_coef;
      // m_se              = Ch.solve(I_p()).diagonal().array().sqrt();
    }
    // } else if (type == 4)
    // {
    // //     UDV.compute((w.asDiagonal() *
    // X).jacobiSvd(ComputeThinU|ComputeThinV));
    // //     MatrixXd             VDi(UDV.matrixV() *
    // //         Dplus(UDV.singularValues().array()).matrix().asDiagonal());
    // //     beta                   = VDi * UDV.matrixU().adjoint() *
    // (z.array() * w.array()).matrix();
    // //     //m_fitted                 = X * m_coef;
    // //     //m_se                     = VDi.rowwise().norm();
    // // } else if (type == 5)
    // // {
    //     eig.compute(XtWX().selfadjointView<Lower>());
    //     MatrixXd   VDi(eig.eigenvectors() *
    //         Dplus(eig.eigenvalues().array()).sqrt().matrix().asDiagonal());
    //     beta         = VDi * VDi.adjoint() * X.adjoint() * (z.array() *
    //     w.array()).matrix();
    //     //m_fitted       = X * m_coef;
    //     //m_se           = VDi.rowwise().norm();
    // }
  }

  virtual void save_se() {
    if (type == 0) {
      if (rank == nvars) {  // full rank case
        se = Pmat * PQR.matrixQR()
                        .topRows(nvars)
                        .triangularView<Upper>()
                        .solve(MatrixXd::Identity(nvars, nvars))
                        .rowwise()
                        .norm();
        return;
      } else {
        // create fitted values from effects
        // (can't use X*m_coef if X is rank-deficient)
        se.head(rank) = Rinv.rowwise().norm();
        se = Pmat * se;
      }
    } else if (type == 1) {
      se = QR.matrixQR()
               .topRows(nvars)
               .triangularView<Upper>()
               .solve(MatrixXd::Identity(nvars, nvars))
               .rowwise()
               .norm();
    } else if (type == 2) {
      se =
          Ch.matrixL().solve(MatrixXd::Identity(nvars, nvars)).colwise().norm();
    } else if (type == 3) {
      se =
          ChD.solve(MatrixXd::Identity(nvars, nvars)).diagonal().array().sqrt();
    } else if (type == 4) {
      if (rank == nvars) {  // full rank case
        se = Pmat * FPQR.matrixQR()
                        .topRows(nvars)
                        .triangularView<Upper>()
                        .solve(MatrixXd::Identity(nvars, nvars))
                        .rowwise()
                        .norm();
        return;
      } else {
        // create fitted values from effects
        // (can't use X*m_coef if X is rank-deficient)
        se.head(rank) = Rinv.rowwise().norm();
        se = Pmat * se;
      }
    } else if (type == 5) {
      Rinv = (bSVD.solve(MatrixXd::Identity(nvars, nvars)));
      se = Rinv.rowwise().norm();
    }
  }

 public:
  glm(Map<MatrixXd> const& X_, Map<VectorXd> const& Y_,
      Map<VectorXd> const& weights_, Map<VectorXd> const& offset_,
      NumericVector (*variance_fun_)(Eigen::VectorXd const&),
      NumericVector (*mu_eta_fun_)(Eigen::VectorXd const&),
      NumericVector (*linkinv_)(Eigen::VectorXd const&),
      NumericVector (*dev_resids_fun_)(Map<VectorXd> const&,
                                       Eigen::VectorXd const&,
                                       Map<VectorXd> const&),
      bool (*valideta_)(Eigen::VectorXd const&),
      bool (*validmu_)(Eigen::VectorXd const&), double tol_ = 1e-6,
      int maxit_ = 100, int type_ = 1, bool is_big_matrix_ = false)
      : GlmBase<Eigen::VectorXd, Eigen::MatrixXd>(X_.rows(), X_.cols(), tol_,
                                                  maxit_),
        X(X_),
        Y(Y_),
        weights(weights_),
        offset(offset_),
        // X(X_.data(), X_.rows(), X_.cols()),
        // Y(Y_.data(), Y_.size()),
        variance_fun(variance_fun_),
        mu_eta_fun(mu_eta_fun_),
        linkinv(linkinv_),
        dev_resids_fun(dev_resids_fun_),
        valideta(valideta_),
        validmu(validmu_),
        tol(tol_),
        maxit(maxit_),
        type(type_),
        is_big_matrix(is_big_matrix_) {}

  // must set params to starting vals
  void init_parms(Map<VectorXd> const& start_, Map<VectorXd> const& mu_,
                  Map<VectorXd> const& eta_) {
    beta = start_;
    eta = eta_;
    mu = mu_;
    // eta.array() += offset.array();

    // update_var_mu();

    // update_mu_eta();

    // update_mu();

    update_dev_resids();

    // update_z();

    // update_w();

    rank = nvars;
  }

  virtual VectorXd get_beta() {
    if (type == 0 || type == 4) {
      if (rank != nvars) {
        // beta.head(rank)                 = Rinv * effects.head(rank);
        // beta = Pmat * beta;
      }
    }

    return beta;
  }

  virtual VectorXd get_weights() { return weights; }
  virtual int get_rank() { return rank; }
};

#endif  // GLM_H
