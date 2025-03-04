#include "fastcpd_class.h"
#include "ref_fastglm_fit_glm.h"
#include "ref_tseries.h"

using ::arma::as_scalar;
using ::arma::colvec;
using ::arma::cube;
using ::arma::eye;
using ::arma::log_det_sympd;
using ::arma::mat;
using ::arma::ones;
using ::arma::reverse;
using ::arma::rowvec;
using ::arma::trace;
using ::arma::ucolvec;
using ::arma::zeros;
using ::Rcpp::as;
using ::Rcpp::Function;
using ::Rcpp::List;
using ::Rcpp::Named;
using ::Rcpp::Nullable;
using ::Rcpp::wrap;
using ::std::pow;
using ::std::string;
using ::std::string_view;
using ::std::unique_ptr;
using ::std::unordered_map;
using ::std::vector;

using ::Rcpp::Environment;
using ::Rcpp::NumericVector;
using ::Rcpp::S4;

namespace fastcpd::classes {

colvec Fastcpd::GetGradientArma(const unsigned int segment_start,
                                const unsigned int segment_end,
                                const colvec& theta) {
  const mat data_segment = data_.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  mat reversed_data = reverse(data_segment, 0);
  colvec reversed_theta = reverse(theta);
  if (segment_length < max(order_) + 1) {
    return ones(theta.n_elem);
  }
  colvec variance_term = zeros(segment_length);
  for (unsigned int i = max(order_); i < segment_length; i++) {
    variance_term(i) = data_segment(i, 0) -
                       dot(reversed_theta.rows(order_(1) + 1, sum(order_)),
                           data_segment.rows(i - order_(0), i - 1)) -
                       dot(reversed_theta.rows(1, order_(1)),
                           variance_term.rows(i - order_(1), i - 1));
  }
  colvec reversed_variance_term = reverse(variance_term);
  mat phi_coefficient = zeros(segment_length, order_(0)),
      psi_coefficient = zeros(segment_length, order_(1));
  for (unsigned int i = max(order_); i < segment_length; i++) {
    phi_coefficient.row(i) =
        -reversed_data
             .rows(segment_length - i, segment_length - i + order_(0) - 1)
             .t() -
        reversed_theta.rows(1, order_(1)).t() *
            phi_coefficient.rows(i - order_(1), i - 1);
  }
  for (unsigned int i = order_(1); i < segment_length; i++) {
    psi_coefficient.row(i) =
        -reversed_variance_term
             .rows(segment_length - i, segment_length - i + order_(1) - 1)
             .t() -
        reversed_theta.rows(1, order_(1)).t() *
            psi_coefficient.rows(i - order_(1), i - 1);
  }
  colvec gradient = zeros(sum(order_) + 1);
  gradient.rows(0, order_(0) - 1) =
      phi_coefficient.row(segment_length - 1).t() *
      variance_term(segment_length - 1) / theta(sum(order_));
  gradient.rows(order_(0), sum(order_) - 1) =
      psi_coefficient.row(segment_length - 1).t() *
      variance_term(segment_length - 1) / theta(sum(order_));
  gradient(sum(order_)) = 1.0 / 2.0 / theta(sum(order_)) -
                          pow(variance_term(segment_length - 1), 2) / 2.0 /
                              pow(theta(sum(order_)), 2);
  return gradient;
}

colvec Fastcpd::GetGradientBinomial(const unsigned int segment_start,
                                    const unsigned int segment_end,
                                    const colvec& theta) {
  const mat data_segment = data_.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double y = new_data(0);
  return -(y - 1 / (1 + exp(-as_scalar(x * theta)))) * x.t();
}

colvec Fastcpd::GetGradientCustom(const unsigned int segment_start,
                                  const unsigned int segment_end,
                                  const colvec& theta) {
  return as<colvec>(
      (*cost_gradient_)(data_.rows(segment_start, segment_end), theta));
}

colvec Fastcpd::GetGradientLm(const unsigned int segment_start,
                              const unsigned int segment_end,
                              const colvec& theta) {
  const mat data_segment = data_.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double y = new_data(0);
  return -(y - as_scalar(x * theta)) * x.t();
}

colvec Fastcpd::GetGradientMa(const unsigned int segment_start,
                              const unsigned int segment_end,
                              const colvec& theta) {
  const mat data_segment = data_.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  const unsigned int q = order_(1);
  mat reversed_data = reverse(data_segment, 0);
  colvec reversed_theta = reverse(theta);
  if (segment_length < q + 1) {
    return ones(theta.n_elem);
  }
  colvec variance_term = zeros(segment_length);
  for (unsigned int i = q; i < segment_length; i++) {
    variance_term(i) =
        data_segment(i, 0) -
        dot(reversed_theta.rows(1, q), variance_term.rows(i - q, i - 1));
  }
  colvec reversed_variance_term = reverse(variance_term);
  mat psi_coefficient = zeros(segment_length, q);
  for (unsigned int i = q; i < segment_length; i++) {
    psi_coefficient.row(i) =
        -reversed_variance_term
             .rows(segment_length - i, segment_length - i + q - 1)
             .t() -
        reversed_theta.rows(1, q).t() * psi_coefficient.rows(i - q, i - 1);
  }
  colvec gradient = zeros(q + 1);
  gradient.rows(0, q - 1) = psi_coefficient.row(segment_length - 1).t() *
                            variance_term(segment_length - 1) / theta(q);
  gradient(q) =
      1.0 / 2.0 / theta(q) -
      pow(variance_term(segment_length - 1), 2) / 2.0 / pow(theta(q), 2);
  return gradient;
}

colvec Fastcpd::GetGradientPoisson(const unsigned int segment_start,
                                   const unsigned int segment_end,
                                   const colvec& theta) {
  const mat data_segment = data_.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double y = new_data(0);
  return -(y - exp(as_scalar(x * theta))) * x.t();
}

mat Fastcpd::GetHessianArma(const unsigned int segment_start,
                            const unsigned int segment_end,
                            const colvec& theta) {
  const mat data_segment = data_.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  // TODO(doccstat): Maybe we can store all these computations
  mat reversed_data = reverse(data_segment, 0);
  colvec reversed_theta = reverse(theta);
  if (segment_length < max(order_) + 1) {
    return eye(theta.n_elem, theta.n_elem);
  }
  colvec variance_term = zeros(segment_length);
  for (unsigned int i = max(order_); i < segment_length; i++) {
    variance_term(i) = data_segment(i, 0) -
                       dot(reversed_theta.rows(order_(1) + 1, sum(order_)),
                           data_segment.rows(i - order_(0), i - 1)) -
                       dot(reversed_theta.rows(1, order_(1)),
                           variance_term.rows(i - order_(1), i - 1));
  }
  colvec reversed_variance_term = reverse(variance_term);
  mat phi_coefficient = zeros(segment_length, order_(0)),
      psi_coefficient = zeros(segment_length, order_(1));
  for (unsigned int i = max(order_); i < segment_length; i++) {
    phi_coefficient.row(i) =
        -reversed_data
             .rows(segment_length - i, segment_length - i + order_(0) - 1)
             .t() -
        reversed_theta.rows(1, order_(1)).t() *
            phi_coefficient.rows(i - order_(1), i - 1);
  }
  for (unsigned int i = order_(1); i < segment_length; i++) {
    psi_coefficient.row(i) =
        -reversed_variance_term
             .rows(segment_length - i, segment_length - i + order_(1) - 1)
             .t() -
        reversed_theta.rows(1, order_(1)).t() *
            psi_coefficient.rows(i - order_(1), i - 1);
  }
  mat reversed_coef_phi = reverse(phi_coefficient, 0),
      reversed_coef_psi = reverse(psi_coefficient, 0);
  cube phi_psi_coefficient = zeros(order_(1), order_(0), segment_length),
       psi_psi_coefficient = zeros(order_(1), order_(1), segment_length);
  for (unsigned int i = order_(1); i < segment_length; i++) {
    mat phi_psi_coefficient_part = zeros(order_(1), order_(0)),
        psi_psi_coefficient_part = zeros(order_(1), order_(1));
    for (unsigned int j = 1; j <= order_(1); j++) {
      phi_psi_coefficient_part +=
          phi_psi_coefficient.slice(i - j) * theta(order_(0) - 1 + j);
    }
    phi_psi_coefficient.slice(i) =
        -reversed_coef_phi.rows(segment_length - i,
                                segment_length - i + order_(1) - 1) -
        phi_psi_coefficient_part;
    for (unsigned int j = 1; j <= order_(1); j++) {
      psi_psi_coefficient_part +=
          psi_psi_coefficient.slice(i - j) * theta(order_(0) - 1 + j);
    }
    psi_psi_coefficient.slice(i) =
        -reversed_coef_psi.rows(segment_length - i,
                                segment_length - i + order_(1) - 1) -
        reversed_coef_psi
            .rows(segment_length - i, segment_length - i + order_(1) - 1)
            .t() -
        psi_psi_coefficient_part;
  }
  mat hessian = zeros(sum(order_) + 1, sum(order_) + 1);
  hessian.submat(0, 0, order_(0) - 1, order_(0) - 1) =
      phi_coefficient.row(segment_length - 1).t() *
      phi_coefficient.row(segment_length - 1) / theta(sum(order_));
  hessian.submat(0, order_(0), order_(0) - 1, sum(order_) - 1) =
      (phi_psi_coefficient.slice(segment_length - 1).t() *
           variance_term(segment_length - 1) +
       phi_coefficient.row(segment_length - 1).t() *
           psi_coefficient.row(segment_length - 1)) /
      theta(sum(order_));
  hessian.submat(order_(0), 0, sum(order_) - 1, order_(0) - 1) =
      hessian.submat(0, order_(0), order_(0) - 1, sum(order_) - 1).t();
  hessian.submat(0, sum(order_), order_(0) - 1, sum(order_)) =
      -phi_coefficient.row(segment_length - 1).t() *
      variance_term(segment_length - 1) / theta(sum(order_)) /
      theta(sum(order_));
  hessian.submat(sum(order_), 0, sum(order_), order_(0) - 1) =
      hessian.submat(0, sum(order_), order_(0) - 1, sum(order_)).t();
  hessian.submat(order_(0), order_(0), sum(order_) - 1, sum(order_) - 1) =
      (psi_coefficient.row(segment_length - 1).t() *
           psi_coefficient.row(segment_length - 1) +
       psi_psi_coefficient.slice(segment_length - 1) *
           variance_term(segment_length - 1)) /
      theta(sum(order_));
  hessian.submat(order_(0), sum(order_), sum(order_) - 1, sum(order_)) =
      -psi_coefficient.row(segment_length - 1).t() *
      variance_term(segment_length - 1) / theta(sum(order_)) /
      theta(sum(order_));
  hessian.submat(sum(order_), order_(0), sum(order_), sum(order_) - 1) =
      hessian.submat(order_(0), sum(order_), sum(order_) - 1, sum(order_)).t();
  hessian(sum(order_), sum(order_)) =
      pow(variance_term(segment_length - 1), 2) / pow(theta(sum(order_)), 3) -
      1.0 / 2.0 / pow(theta(sum(order_)), 2);
  return hessian;
}

mat Fastcpd::GetHessianBinomial(const unsigned int segment_start,
                                const unsigned int segment_end,
                                const colvec& theta) {
  const mat data_segment = data_.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double prob = 1 / (1 + exp(-dot(x, theta)));
  return (x.t() * x) * ((1 - prob) * prob);
}

mat Fastcpd::GetHessianCustom(const unsigned int segment_start,
                              const unsigned int segment_end,
                              const colvec& theta) {
  return as<mat>(
      (*cost_hessian_)(data_.rows(segment_start, segment_end), theta));
}

mat Fastcpd::GetHessianLm(const unsigned int segment_start,
                          const unsigned int segment_end, const colvec& theta) {
  const mat data_segment = data_.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  return x.t() * x;
}

mat Fastcpd::GetHessianMa(const unsigned int segment_start,
                          const unsigned int segment_end, const colvec& theta) {
  const mat data_segment = data_.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  const unsigned int q = order_(1);
  // TODO(doccstat): Maybe we can store all these computations
  mat reversed_data = reverse(data_segment, 0);
  colvec reversed_theta = reverse(theta);
  if (segment_length < q + 1) {
    return eye(theta.n_elem, theta.n_elem);
  }
  colvec variance_term = zeros(segment_length);
  for (unsigned int i = q; i < segment_length; i++) {
    variance_term(i) =
        data_segment(i, 0) -
        dot(reversed_theta.rows(1, q), variance_term.rows(i - q, i - 1));
  }
  colvec reversed_variance_term = reverse(variance_term);
  mat psi_coefficient = zeros(segment_length, q);
  for (unsigned int i = q; i < segment_length; i++) {
    psi_coefficient.row(i) =
        -reversed_variance_term
             .rows(segment_length - i, segment_length - i + q - 1)
             .t() -
        reversed_theta.rows(1, q).t() * psi_coefficient.rows(i - q, i - 1);
  }
  mat reversed_coef_psi = reverse(psi_coefficient, 0);
  cube psi_psi_coefficient = zeros(q, q, segment_length);
  for (unsigned int i = q; i < segment_length; i++) {
    mat psi_psi_coefficient_part = zeros(q, q);
    for (unsigned int j = 1; j <= q; j++) {
      psi_psi_coefficient_part +=
          psi_psi_coefficient.slice(i - j) * theta(j - 1);
    }
    psi_psi_coefficient.slice(i) =
        -reversed_coef_psi.rows(segment_length - i,
                                segment_length - i + q - 1) -
        reversed_coef_psi.rows(segment_length - i, segment_length - i + q - 1)
            .t() -
        psi_psi_coefficient_part;
  }
  mat hessian = zeros(q + 1, q + 1);
  hessian.submat(0, 0, q - 1, q - 1) =
      (psi_coefficient.row(segment_length - 1).t() *
           psi_coefficient.row(segment_length - 1) +
       psi_psi_coefficient.slice(segment_length - 1) *
           variance_term(segment_length - 1)) /
      theta(q);
  hessian.submat(0, q, q - 1, q) =
      -psi_coefficient.row(segment_length - 1).t() *
      variance_term(segment_length - 1) / theta(q) / theta(q);
  hessian.submat(q, 0, q, q - 1) = hessian.submat(0, q, q - 1, q).t();
  hessian(q, q) = pow(variance_term(segment_length - 1), 2) / pow(theta(q), 3) -
                  1.0 / 2.0 / pow(theta(q), 2);
  return hessian;
}

mat Fastcpd::GetHessianPoisson(const unsigned int segment_start,
                               const unsigned int segment_end,
                               const colvec& theta) {
  const mat data_segment = data_.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double prob = exp(as_scalar(x * theta));
  // Prevent numerical issues if `prob` is too large.
  return (x.t() * x) * std::min(as_scalar(prob), 1e10);
}

CostResult Fastcpd::GetNllPeltArma(const unsigned int segment_start,
                                   const unsigned int segment_end,
                                   const bool cv,
                                   const Nullable<colvec>& start) {
  const mat data_segment = data_.rows(segment_start, segment_end);
  Environment stats = Environment::namespace_env("stats");
  Function arima = stats["arima"];
  List out =
      arima(Named("x") = data_segment.col(0),
            Named("order") = NumericVector::create(order_(0), 0, order_(1)),
            Named("include.mean") = false);
  colvec par = zeros(sum(order_) + 1);
  par.rows(0, sum(order_) - 1) = as<colvec>(out["coef"]);
  par(sum(order_)) = as<double>(out["sigma2"]);

  return {{par}, {as<colvec>(out["residuals"])}, -as<double>(out["loglik"])};
}

CostResult Fastcpd::GetNllPeltCustom(const unsigned int segment_start,
                                     const unsigned int segment_end,
                                     const bool cv,
                                     const Nullable<colvec>& start) {
  if (cost_gradient_ || cost_hessian_) {
    return GetOptimizedCostResult(segment_start, segment_end);
  } else {
    return {
        {colvec()},
        {colvec()},
        as<double>((*cost_function_)(data_.rows(segment_start, segment_end)))};
  }
}

CostResult Fastcpd::GetNllPeltGarch(const unsigned int segment_start,
                                    const unsigned int segment_end,
                                    const bool cv,
                                    const Nullable<colvec>& start) {
  colvec series = data_.rows(segment_start, segment_end).col(0);
  List out = garch(series, order_);
  return {{as<colvec>(out["coef"])},
          {as<colvec>(out["residuals"])},
          as<double>(out["n.likeli"])};
}

CostResult Fastcpd::GetNllPeltGlm(const unsigned int segment_start,
                                  const unsigned int segment_end, const bool cv,
                                  const Nullable<colvec>& start) {
  const mat data_segment = data_.rows(segment_start, segment_end);
  colvec y = data_segment.col(0);
  List out;
  if (start.isNull()) {
    mat x = data_segment.cols(1, data_segment.n_cols - 1);
    out = fastglm(x, y, family_);
  } else {
    mat x = data_segment.cols(1, data_segment.n_cols - 1);
    out = fastglm(x, y, family_, start);
  }
  colvec par = as<colvec>(out["coefficients"]);
  colvec residuals = as<colvec>(out["residuals"]);
  double value = out["deviance"];
  return {{par}, {residuals}, value / 2};
}

CostResult Fastcpd::GetNllPeltLasso(const unsigned int segment_start,
                                    const unsigned int segment_end,
                                    const bool cv,
                                    const Nullable<colvec>& start) {
  if (segment_start == segment_end) {
    return {{zeros(data_.n_cols - 1)}, {zeros(1)}, 0};
  }
  if (cv) {
    const mat data_segment = data_.rows(segment_start, segment_end);
    Environment glmnet = Environment::namespace_env("glmnet"),
                stats = Environment::namespace_env("stats");
    Function cv_glmnet = glmnet["cv.glmnet"],
             predict_glmnet = glmnet["predict.glmnet"],
             deviance = stats["deviance"];
    List out = cv_glmnet(data_segment.cols(1, data_segment.n_cols - 1),
                         data_segment.col(0), Named("family") = "gaussian");
    colvec index_vec = as<colvec>(out["index"]),
           values = as<colvec>(deviance(out["glmnet.fit"]));
    S4 out_coef =
        predict_glmnet(out["glmnet.fit"], Named("s") = out["lambda.1se"],
                       Named("type") = "coefficients", Named("exact") = false);
    colvec glmnet_i = as<colvec>(out_coef.slot("i"));
    colvec glmnet_x = as<colvec>(out_coef.slot("x"));
    colvec par = zeros(data_segment.n_cols - 1);
    for (unsigned int i = 1; i < glmnet_i.n_elem; i++) {
      par(glmnet_i(i) - 1) = glmnet_x(i);
    }
    return {{par}, {mat()}, values(index_vec(1) - 1)};
  } else {
    const mat data_segment = data_.rows(segment_start, segment_end);
    Environment stats = Environment::namespace_env("stats"),
                glmnet = Environment::namespace_env("glmnet");
    Function deviance = stats["deviance"], glmnet_ = glmnet["glmnet"],
             predict_glmnet = glmnet["predict.glmnet"];
    List out = glmnet_(data_segment.cols(1, data_segment.n_cols - 1),
                       data_segment.col(0), Named("family") = "gaussian",
                       Named("lambda") = lasso_penalty_base_ /
                                         sqrt(segment_end - segment_start + 1));
    S4 out_par = out["beta"];
    colvec par_i = as<colvec>(out_par.slot("i"));
    colvec par_x = as<colvec>(out_par.slot("x"));
    colvec par = zeros(data_segment.n_cols - 1);
    for (unsigned int i = 0; i < par_i.n_elem; i++) {
      par(par_i(i)) = par_x(i);
    }
    double value = as<double>(deviance(out));
    colvec fitted_values = as<colvec>(predict_glmnet(
        out, data_segment.cols(1, data_segment.n_cols - 1),
        Named("s") =
            lasso_penalty_base_ / sqrt(segment_end - segment_start + 1)));
    colvec residuals = data_segment.col(0) - fitted_values;
    return {{par}, {residuals}, value / 2};
  }
}

CostResult Fastcpd::GetNllPeltMean(const unsigned int segment_start,
                                   const unsigned int segment_end,
                                   const bool cv,
                                   const Nullable<colvec>& start) {
  double two_norm = 0;
  for (unsigned int i = 0; i < parameters_count_; i++) {
    two_norm +=
        (data_c_.at(segment_end + 1, i) - data_c_.at(segment_start, i)) *
        (data_c_.at(segment_end + 1, i) - data_c_.at(segment_start, i));
  }
  const unsigned int segment_length = segment_end - segment_start + 1;
  return {{zeros<colvec>(parameters_count_)},
          {zeros<mat>(segment_length, parameters_count_)},
          ((data_c_.at(segment_end + 1, parameters_count_) -
            data_c_.at(segment_start, parameters_count_)) -
           two_norm / segment_length) /
              2.0};
}

CostResult Fastcpd::GetNllPeltMeanVariance(const unsigned int segment_start,
                                           const unsigned int segment_end,
                                           const bool cv,
                                           const Nullable<colvec>& start) {
  rowvec data_diff = data_c_.row(segment_end + 1) - data_c_.row(segment_start);
  const unsigned int segment_length = segment_end - segment_start + 1;

  double det_value =
      det((reshape(data_diff.subvec(data_n_dims_, parameters_count_ - 1),
                   data_n_dims_, data_n_dims_) -
           (data_diff.subvec(0, data_n_dims_ - 1)).t() *
               (data_diff.subvec(0, data_n_dims_ - 1)) / segment_length) /
          segment_length);
  if (segment_length <= data_n_dims_) {
    unsigned int approximate_segment_start;
    unsigned int approximate_segment_end;
    if (segment_start >= data_n_dims_) {
      approximate_segment_start = segment_start - data_n_dims_;
    } else {
      approximate_segment_start = 0;
    }
    if (segment_end < data_n_rows_ - data_n_dims_) {
      approximate_segment_end = segment_end + data_n_dims_;
    } else {
      approximate_segment_end = data_n_rows_ - 1;
    }
    data_diff = data_c_.row(approximate_segment_end + 1) -
                data_c_.row(approximate_segment_start);
    det_value =
        det((reshape(data_diff.subvec(data_n_dims_, parameters_count_ - 1),
                     data_n_dims_, data_n_dims_) -
             (data_diff.subvec(0, data_n_dims_ - 1)).t() *
                 (data_diff.subvec(0, data_n_dims_ - 1)) /
                 (approximate_segment_end - approximate_segment_start + 1)) /
            (approximate_segment_end - approximate_segment_start + 1));
  }

  return {{zeros<colvec>(parameters_count_)},
          {mat()},
          log(det_value) * segment_length / 2.0};
}

CostResult Fastcpd::GetNllPeltMgaussian(const unsigned int segment_start,
                                        const unsigned int segment_end,
                                        const bool cv,
                                        const Nullable<colvec>& start) {
  const mat data_segment = data_.rows(segment_start, segment_end);
  mat x =
      data_segment.cols(regression_response_count_, data_segment.n_cols - 1);
  mat y = data_segment.cols(0, regression_response_count_ - 1);
  mat x_t_x;

  if (data_segment.n_rows <=
      data_segment.n_cols - regression_response_count_ + 1) {
    x_t_x = eye<mat>(data_segment.n_cols - regression_response_count_,
                     data_segment.n_cols - regression_response_count_);
  } else {
    x_t_x = x.t() * x;
  }

  mat par = solve(x_t_x, x.t()) * y;
  mat residuals = y - x * par;
  double value = regression_response_count_ * std::log(2.0 * M_PI) +
                 log_det_sympd(variance_estimate_);
  value *= data_segment.n_rows;
  value += trace(solve(variance_estimate_, residuals.t() * residuals));
  return {{par}, {residuals}, value / 2};
}

CostResult Fastcpd::GetNllPeltVariance(const unsigned int segment_start,
                                       const unsigned int segment_end,
                                       const bool cv,
                                       const Nullable<colvec>& start) {
  const unsigned int segment_length = segment_end - segment_start + 1;

  double det_value = det(
      arma::reshape(data_c_.row(segment_end + 1) - data_c_.row(segment_start),
                    data_n_dims_, data_n_dims_) /
      segment_length);
  if (segment_length < data_n_dims_) {
    unsigned int approximate_segment_start;
    unsigned int approximate_segment_end;
    if (segment_start >= data_n_dims_) {
      approximate_segment_start = segment_start - data_n_dims_;
    } else {
      approximate_segment_start = 0;
    }
    if (segment_end < data_n_rows_ - data_n_dims_) {
      approximate_segment_end = segment_end + data_n_dims_;
    } else {
      approximate_segment_end = data_n_rows_ - 1;
    }
    det_value = det(arma::reshape(data_c_.row(approximate_segment_end + 1) -
                                      data_c_.row(approximate_segment_start),
                                  data_n_dims_, data_n_dims_) /
                    (approximate_segment_end - approximate_segment_start + 1));
  }

  return {{zeros<mat>(data_n_dims_, data_n_dims_)},
          {mat()},
          log(det_value) * segment_length / 2.0};
}

double Fastcpd::GetNllSenArma(const unsigned int segment_start,
                              const unsigned int segment_end, colvec theta) {
  mat data_segment = data_.rows(segment_start, segment_end);
  colvec reversed_theta = reverse(theta);
  if (data_segment.n_rows < max(order_) + 1) {
    return 0;
  }
  colvec variance_term = zeros(data_segment.n_rows);
  for (unsigned int i = max(order_); i < data_segment.n_rows; i++) {
    variance_term(i) = data_segment(i, 0) -
                       dot(reversed_theta.rows(order_(1) + 1, sum(order_)),
                           data_segment.rows(i - order_(0), i - 1)) -
                       dot(reversed_theta.rows(1, order_(1)),
                           variance_term.rows(i - order_(1), i - 1));
  }
  return (std::log(2.0 * M_PI) + std::log(theta(sum(order_)))) *
             (data_segment.n_rows - 2) / 2.0 +
         dot(variance_term, variance_term) / 2.0 / theta(sum(order_));
}

double Fastcpd::GetNllSenBinomial(const unsigned int segment_start,
                                  const unsigned int segment_end,
                                  colvec theta) {
  mat data_segment = data_.rows(segment_start, segment_end);
  colvec y = data_segment.col(0);
  // Calculate negative log likelihood in binomial family
  mat x = data_segment.cols(1, data_segment.n_cols - 1);
  colvec u = x * theta;
  return accu(-y % u + arma::log(1 + exp(u)));
}

double Fastcpd::GetNllSenCustom(const unsigned int segment_start,
                                const unsigned int segment_end, colvec theta) {
  return as<double>(
      (*cost_function_)(data_.rows(segment_start, segment_end), theta));
}

double Fastcpd::GetNllSenLasso(const unsigned int segment_start,
                               const unsigned int segment_end, colvec theta) {
  mat data_segment = data_.rows(segment_start, segment_end);
  colvec y = data_segment.col(0);
  mat x = data_segment.cols(1, data_segment.n_cols - 1);
  return accu(square(y - x * theta)) / 2 +
         lasso_penalty_base_ / sqrt(segment_end - segment_start + 1) *
             accu(abs(theta));
}

double Fastcpd::GetNllSenLm(const unsigned int segment_start,
                            const unsigned int segment_end, colvec theta) {
  mat data_segment = data_.rows(segment_start, segment_end);
  colvec y = data_segment.col(0);
  mat x = data_segment.cols(1, data_segment.n_cols - 1);
  return accu(square(y - x * theta)) / 2;
}

double Fastcpd::GetNllSenMa(const unsigned int segment_start,
                            const unsigned int segment_end, colvec theta) {
  mat data_segment = data_.rows(segment_start, segment_end);
  const unsigned int q = order_(1);
  colvec reversed_theta = reverse(theta);
  if (data_segment.n_rows < q + 1) {
    return 0;
  }
  colvec variance_term = zeros(data_segment.n_rows);
  for (unsigned int i = q; i < data_segment.n_rows; i++) {
    variance_term(i) =
        data_segment(i, 0) -
        dot(reversed_theta.rows(1, q), variance_term.rows(i - q, i - 1));
  }
  return (std::log(2.0 * M_PI) + std::log(theta(q))) *
             (data_segment.n_rows - 2) / 2.0 +
         dot(variance_term, variance_term) / 2.0 / theta(q);
}

double Fastcpd::GetNllSenPoisson(const unsigned int segment_start,
                                 const unsigned int segment_end, colvec theta) {
  mat data_segment = data_.rows(segment_start, segment_end);
  colvec y = data_segment.col(0);
  mat x = data_segment.cols(1, data_segment.n_cols - 1);
  colvec u = x * theta;
  colvec y_factorial(y.n_elem);
  for (unsigned int i = 0; i < y.n_elem; i++) {
    double log_factorial = 0;
    for (int j = 1; j <= y(i); ++j) {
      log_factorial += std::log(j);
    }
    y_factorial(i) = log_factorial;
  }
  return accu(-y % u + exp(u) + y_factorial);
}

}  // namespace fastcpd::classes
