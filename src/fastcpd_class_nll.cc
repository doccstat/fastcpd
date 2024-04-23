#include "fastcpd_classes.h"

namespace fastcpd::classes {

colvec Fastcpd::get_gradient_arma(
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) {
  const mat data_segment = data.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  mat reversed_data = reverse(data_segment, 0);
  colvec reversed_theta = reverse(theta);
  if (segment_length < max(order) + 1) {
    return ones(theta.n_elem);
  }
  colvec variance_term = zeros(segment_length);
  for (unsigned int i = max(order); i < segment_length; i++) {
    variance_term(i) = data_segment(i, 0) - dot(
        reversed_theta.rows(order(1) + 1, sum(order)),
        data_segment.rows(i - order(0), i - 1)
      ) - dot(
        reversed_theta.rows(1, order(1)),
        variance_term.rows(i - order(1), i - 1)
      );
  }
  colvec reversed_variance_term = reverse(variance_term);
  mat phi_coefficient = zeros(segment_length, order(0)),
      psi_coefficient = zeros(segment_length, order(1));
  for (unsigned int i = max(order); i < segment_length; i++) {
    phi_coefficient.row(i) = -reversed_data.rows(
      segment_length - i, segment_length - i + order(0) - 1
    ).t() - reversed_theta.rows(1, order(1)).t() *
    phi_coefficient.rows(i - order(1), i - 1);
  }
  for (unsigned int i = order(1); i < segment_length; i++) {
    psi_coefficient.row(i) = -reversed_variance_term.rows(
        segment_length - i, segment_length - i + order(1) - 1
      ).t() - reversed_theta.rows(1, order(1)).t() *
      psi_coefficient.rows(i - order(1), i - 1);
  }
  colvec gradient = zeros(sum(order) + 1);
  gradient.rows(0, order(0) - 1) =
    phi_coefficient.row(segment_length - 1).t() *
    variance_term(segment_length - 1) / theta(sum(order));
  gradient.rows(order(0), sum(order) - 1) =
    psi_coefficient.row(segment_length - 1).t() *
    variance_term(segment_length - 1) / theta(sum(order));
  gradient(sum(order)) = 1.0 / 2.0 / theta(sum(order)) -
    std::pow(variance_term(segment_length - 1), 2) / 2.0 /
    std::pow(theta(sum(order)), 2);
  return gradient;
}

colvec Fastcpd::get_gradient_binomial(
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) {
  const mat data_segment = data.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double y = new_data(0);
  return - (y - 1 / (1 + exp(-as_scalar(x * theta)))) * x.t();
}

colvec Fastcpd::get_gradient_lm(
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) {
  const mat data_segment = data.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double y = new_data(0);
  return - (y - as_scalar(x * theta)) * x.t();
}

colvec Fastcpd::get_gradient_ma(
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) {
  const mat data_segment = data.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  const unsigned int q = order(1);
  mat reversed_data = reverse(data_segment, 0);
  colvec reversed_theta = reverse(theta);
  if (segment_length < q + 1) {
    return ones(theta.n_elem);
  }
  colvec variance_term = zeros(segment_length);
  for (unsigned int i = q; i < segment_length; i++) {
    variance_term(i) = data_segment(i, 0) -
      dot(reversed_theta.rows(1, q), variance_term.rows(i - q, i - 1));
  }
  colvec reversed_variance_term = reverse(variance_term);
  mat psi_coefficient = zeros(segment_length, q);
  for (unsigned int i = q; i < segment_length; i++) {
    psi_coefficient.row(i) = -reversed_variance_term.rows(
        segment_length - i, segment_length - i + q - 1
      ).t() - reversed_theta.rows(1, q).t() *
      psi_coefficient.rows(i - q, i - 1);
  }
  colvec gradient = zeros(q + 1);
  gradient.rows(0, q - 1) =
    psi_coefficient.row(segment_length - 1).t() *
    variance_term(segment_length - 1) / theta(q);
  gradient(q) = 1.0 / 2.0 / theta(q) -
    std::pow(variance_term(segment_length - 1), 2) / 2.0 /
    std::pow(theta(q), 2);
  return gradient;
}

colvec Fastcpd::get_gradient_poisson(
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) {
  const mat data_segment = data.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double y = new_data(0);
  return - (y - exp(as_scalar(x * theta))) * x.t();
}

mat Fastcpd::get_hessian_arma(
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) {
  const mat data_segment = data.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  // TODO(doccstat): Maybe we can store all these computations
  mat reversed_data = reverse(data_segment, 0);
  colvec reversed_theta = reverse(theta);
  if (segment_length < max(order) + 1) {
    return eye(theta.n_elem, theta.n_elem);
  }
  colvec variance_term = zeros(segment_length);
  for (unsigned int i = max(order); i < segment_length; i++) {
    variance_term(i) = data_segment(i, 0) - dot(
        reversed_theta.rows(order(1) + 1, sum(order)),
        data_segment.rows(i - order(0), i - 1)
      ) - dot(
        reversed_theta.rows(1, order(1)),
        variance_term.rows(i - order(1), i - 1)
      );
  }
  colvec reversed_variance_term = reverse(variance_term);
  mat phi_coefficient = zeros(segment_length, order(0)),
      psi_coefficient = zeros(segment_length, order(1));
  for (unsigned int i = max(order); i < segment_length; i++) {
    phi_coefficient.row(i) = -reversed_data.rows(
      segment_length - i, segment_length - i + order(0) - 1
    ).t() - reversed_theta.rows(1, order(1)).t() *
    phi_coefficient.rows(i - order(1), i - 1);
  }
  for (unsigned int i = order(1); i < segment_length; i++) {
    psi_coefficient.row(i) = -reversed_variance_term.rows(
        segment_length - i, segment_length - i + order(1) - 1
      ).t() - reversed_theta.rows(1, order(1)).t() *
      psi_coefficient.rows(i - order(1), i - 1);
  }
  mat reversed_coef_phi = reverse(phi_coefficient, 0),
      reversed_coef_psi = reverse(psi_coefficient, 0);
  cube phi_psi_coefficient = zeros(order(1), order(0), segment_length),
        psi_psi_coefficient = zeros(order(1), order(1), segment_length);
  for (unsigned int i = order(1); i < segment_length; i++) {
    mat phi_psi_coefficient_part = zeros(order(1), order(0)),
        psi_psi_coefficient_part = zeros(order(1), order(1));
    for (unsigned int j = 1; j <= order(1); j++) {
      phi_psi_coefficient_part +=
        phi_psi_coefficient.slice(i - j) * theta(order(0) - 1 + j);
    }
    phi_psi_coefficient.slice(i) = -reversed_coef_phi.rows(
      segment_length - i, segment_length - i + order(1) - 1
    ) - phi_psi_coefficient_part;
    for (unsigned int j = 1; j <= order(1); j++) {
      psi_psi_coefficient_part +=
        psi_psi_coefficient.slice(i - j) * theta(order(0) - 1 + j);
    }
    psi_psi_coefficient.slice(i) = -reversed_coef_psi.rows(
        segment_length - i, segment_length - i + order(1) - 1
      ) - reversed_coef_psi.rows(
        segment_length - i, segment_length - i + order(1) - 1
      ).t() - psi_psi_coefficient_part;
  }
  mat hessian = zeros(sum(order) + 1, sum(order) + 1);
  hessian.submat(0, 0, order(0) - 1, order(0) - 1) =
    phi_coefficient.row(segment_length - 1).t() *
    phi_coefficient.row(segment_length - 1) / theta(sum(order));
  hessian.submat(0, order(0), order(0) - 1, sum(order) - 1) = (
    phi_psi_coefficient.slice(segment_length - 1).t() *
      variance_term(segment_length - 1) +
      phi_coefficient.row(segment_length - 1).t() *
      psi_coefficient.row(segment_length - 1)
  ) / theta(sum(order));
  hessian.submat(order(0), 0, sum(order) - 1, order(0) - 1) =
    hessian.submat(0, order(0), order(0) - 1, sum(order) - 1).t();
  hessian.submat(0, sum(order), order(0) - 1, sum(order)) =
    -phi_coefficient.row(segment_length - 1).t() *
    variance_term(segment_length - 1) / theta(sum(order)) / theta(sum(order));
  hessian.submat(sum(order), 0, sum(order), order(0) - 1) =
    hessian.submat(0, sum(order), order(0) - 1, sum(order)).t();
  hessian.submat(order(0), order(0), sum(order) - 1, sum(order) - 1) = (
    psi_coefficient.row(segment_length - 1).t() *
    psi_coefficient.row(segment_length - 1) +
    psi_psi_coefficient.slice(segment_length - 1) *
    variance_term(segment_length - 1)
  ) / theta(sum(order));
  hessian.submat(order(0), sum(order), sum(order) - 1, sum(order)) =
    -psi_coefficient.row(segment_length - 1).t() *
    variance_term(segment_length - 1) / theta(sum(order)) / theta(sum(order));
  hessian.submat(sum(order), order(0), sum(order), sum(order) - 1) =
    hessian.submat(order(0), sum(order), sum(order) - 1, sum(order)).t();
  hessian(sum(order), sum(order)) =
    std::pow(variance_term(segment_length - 1), 2) /
    std::pow(theta(sum(order)), 3) -
    1.0 / 2.0 / std::pow(theta(sum(order)), 2);
  return hessian;
}

mat Fastcpd::get_hessian_binomial(
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) {
  const mat data_segment = data.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double prob = 1 / (1 + exp(-dot(x, theta)));
  return (x.t() * x) * ((1 - prob) * prob);
}

mat Fastcpd::get_hessian_lm(
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) {
  const mat data_segment = data.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  return x.t() * x;
}

mat Fastcpd::get_hessian_ma(
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) {
  const mat data_segment = data.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  const unsigned int q = order(1);
  // TODO(doccstat): Maybe we can store all these computations
  mat reversed_data = reverse(data_segment, 0);
  colvec reversed_theta = reverse(theta);
  if (segment_length < q + 1) {
    return eye(theta.n_elem, theta.n_elem);
  }
  colvec variance_term = zeros(segment_length);
  for (unsigned int i = q; i < segment_length; i++) {
    variance_term(i) = data_segment(i, 0) - dot(
      reversed_theta.rows(1, q), variance_term.rows(i - q, i - 1)
    );
  }
  colvec reversed_variance_term = reverse(variance_term);
  mat psi_coefficient = zeros(segment_length, q);
  for (unsigned int i = q; i < segment_length; i++) {
    psi_coefficient.row(i) = -reversed_variance_term.rows(
        segment_length - i, segment_length - i + q - 1
      ).t() - reversed_theta.rows(1, q).t() *
      psi_coefficient.rows(i - q, i - 1);
  }
  mat reversed_coef_psi = reverse(psi_coefficient, 0);
  cube psi_psi_coefficient = zeros(q, q, segment_length);
  for (unsigned int i = q; i < segment_length; i++) {
    mat psi_psi_coefficient_part = zeros(q, q);
    for (unsigned int j = 1; j <= q; j++) {
      psi_psi_coefficient_part +=
        psi_psi_coefficient.slice(i - j) * theta(j - 1);
    }
    psi_psi_coefficient.slice(i) = -reversed_coef_psi.rows(
        segment_length - i, segment_length - i + q - 1
      ) - reversed_coef_psi.rows(
        segment_length - i, segment_length - i + q - 1
      ).t() - psi_psi_coefficient_part;
  }
  mat hessian = zeros(q + 1, q + 1);
  hessian.submat(0, 0, q - 1, q - 1) = (
    psi_coefficient.row(segment_length - 1).t() *
    psi_coefficient.row(segment_length - 1) +
    psi_psi_coefficient.slice(segment_length - 1) *
    variance_term(segment_length - 1)
  ) / theta(q);
  hessian.submat(0, q, q - 1, q) =
    -psi_coefficient.row(segment_length - 1).t() *
    variance_term(segment_length - 1) / theta(q) / theta(q);
  hessian.submat(q, 0, q, q - 1) = hessian.submat(0, q, q - 1, q).t();
  hessian(q, q) =
    std::pow(variance_term(segment_length - 1), 2) /
    std::pow(theta(q), 3) -
    1.0 / 2.0 / std::pow(theta(q), 2);
  return hessian;
}

mat Fastcpd::get_hessian_poisson(
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) {
  const mat data_segment = data.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double prob = exp(as_scalar(x * theta));
  // Prevent numerical issues if `prob` is too large.
  return (x.t() * x) * std::min(as_scalar(prob), 1e10);
}

CostResult Fastcpd::get_nll_arma(
  const unsigned int segment_start,
  const unsigned int segment_end
) {
  const mat data_segment = data.rows(segment_start, segment_end);
  Environment stats = Environment::namespace_env("stats");
  Function arima = stats["arima"];
  List out = arima(
    Named("x") = data_segment.col(0),
    Named("order") = NumericVector::create(order(0), 0, order(1)),
    Named("include.mean") = false
  );
  colvec par = zeros(sum(order) + 1);
  par.rows(0, sum(order) - 1) = as<colvec>(out["coef"]);
  par(sum(order)) = as<double>(out["sigma2"]);

  return {{par}, {as<vec>(out["residuals"])}, -as<double>(out["loglik"])};
}

CostResult Fastcpd::get_nll_glm(
  const unsigned int segment_start,
  const unsigned int segment_end,
  Nullable<colvec> start
) {
  const mat data_segment = data.rows(segment_start, segment_end);
  vec y = data_segment.col(0);
  Environment fastglm = Environment::namespace_env("fastglm");
  Function fastglm_ = fastglm["fastglm"];
  List out;
  if (start.isNull()) {
    mat x = data_segment.cols(1, data_segment.n_cols - 1);
    out = fastglm_(x, y, family);
  } else {
    colvec start_ = as<colvec>(start);
    mat x = data_segment.cols(1, data_segment.n_cols - 1);
    out = fastglm_(x, y, family, Named("start") = start_);
  }
  vec par = as<vec>(out["coefficients"]);
  vec residuals = as<vec>(out["residuals"]);
  double value = out["deviance"];
  return {{par}, {residuals}, value / 2};
}

CostResult Fastcpd::get_nll_lasso_cv(
  const unsigned int segment_start,
  const unsigned int segment_end
) {
  const mat data_segment = data.rows(segment_start, segment_end);
  Environment glmnet = Environment::namespace_env("glmnet"),
               stats = Environment::namespace_env("stats");
  Function cv_glmnet = glmnet["cv.glmnet"],
      predict_glmnet = glmnet["predict.glmnet"],
            deviance = stats["deviance"];
  List out = cv_glmnet(
    data_segment.cols(1, data_segment.n_cols - 1),
    data_segment.col(0),
    Named("family") = "gaussian"
  );
  colvec index_vec = as<colvec>(out["index"]),
            values = as<colvec>(deviance(out["glmnet.fit"]));
  S4 out_coef = predict_glmnet(
    out["glmnet.fit"],
    Named("s") = out["lambda.1se"],
    Named("type") = "coefficients",
    Named("exact") = false
  );
  vec glmnet_i = as<vec>(out_coef.slot("i"));
  vec glmnet_x = as<vec>(out_coef.slot("x"));
  vec par = zeros(data_segment.n_cols - 1);
  for (unsigned int i = 1; i < glmnet_i.n_elem; i++) {
    par(glmnet_i(i) - 1) = glmnet_x(i);
  }
  return {{par}, {mat()}, values(index_vec(1) - 1)};
}

CostResult Fastcpd::get_nll_lasso_wo_cv(
  const unsigned int segment_start,
  const unsigned int segment_end,
  const double lambda
) {
  const mat data_segment = data.rows(segment_start, segment_end);
  Environment stats = Environment::namespace_env("stats"),
             glmnet = Environment::namespace_env("glmnet");
  Function deviance = stats["deviance"], glmnet_ = glmnet["glmnet"],
     predict_glmnet = glmnet["predict.glmnet"];
  List out = glmnet_(
    data_segment.cols(1, data_segment.n_cols - 1), data_segment.col(0),
    Named("family") = "gaussian", Named("lambda") = lambda
  );
  S4 out_par = out["beta"];
  vec par_i = as<vec>(out_par.slot("i"));
  vec par_x = as<vec>(out_par.slot("x"));
  vec par = zeros(data_segment.n_cols - 1);
  for (unsigned int i = 0; i < par_i.n_elem; i++) {
    par(par_i(i)) = par_x(i);
  }
  double value = as<double>(deviance(out));
  vec fitted_values = as<vec>(
    predict_glmnet(
      out, data_segment.cols(1, data_segment.n_cols - 1), Named("s") = lambda
    )
  );
  vec residuals = data_segment.col(0) - fitted_values;
  return {{par}, {residuals}, value / 2};
}

CostResult Fastcpd::get_nll_mean(
  const unsigned int segment_start,
  const unsigned int segment_end
) {
  const rowvec data_diff =
    zero_data.row(segment_end + 1) - zero_data.row(segment_start);
  const unsigned int segment_length = segment_end - segment_start + 1;
  return {
    {zeros<colvec>(p)},  // # nocov
    {zeros<mat>(segment_length, p)},
    std::log(2.0 * M_PI) * zero_data.n_cols + log_det_sympd(variance_estimate) *
      (segment_length) / 2.0 + (
      data_diff(p) - dot(
        data_diff.subvec(0, p - 1), data_diff.subvec(0, p - 1)
      ) / segment_length
    ) / 2.0
  };
}

CostResult Fastcpd::get_nll_meanvariance(
  const unsigned int segment_start,
  const unsigned int segment_end
) {
  rowvec data_diff =
    zero_data.row(segment_end + 1) - zero_data.row(segment_start);
  const unsigned int segment_length = segment_end - segment_start + 1;

  double det_value = det((
    reshape(data_diff.subvec(d, p - 1), d, d) - (
      data_diff.subvec(0, d - 1)).t() * (data_diff.subvec(0, d - 1)
    ) / segment_length
  ) / segment_length);
  if (segment_length <= d) {
    unsigned int approximate_segment_start;
    unsigned int approximate_segment_end;
    if (segment_start >= d) {
      approximate_segment_start = segment_start - d;
    } else {
      approximate_segment_start = 0;
    }
    if (segment_end < data_n_rows - d) {
      approximate_segment_end = segment_end + d;
    } else {
      approximate_segment_end = data_n_rows - 1;
    }
    data_diff = zero_data.row(approximate_segment_end + 1) -
      zero_data.row(approximate_segment_start);
    det_value = det((
    reshape(data_diff.subvec(d, p - 1), d, d) - (
      data_diff.subvec(0, d - 1)).t() * (data_diff.subvec(0, d - 1)
    ) / (approximate_segment_end - approximate_segment_start + 1)
  ) / (approximate_segment_end - approximate_segment_start + 1));
  }

  return {
    {zeros<colvec>(p)},
    {mat()},
    (d * std::log(2.0 * M_PI) + d + log(det_value)) * (segment_length) / 2.0
  };
}

CostResult Fastcpd::get_nll_mgaussian(
  const unsigned int segment_start,
  const unsigned int segment_end
) {
  const mat data_segment = data.rows(segment_start, segment_end);
  mat x = data_segment.cols(p_response, data_segment.n_cols - 1);
  mat y = data_segment.cols(0, p_response - 1);
  mat x_t_x;

  if (data_segment.n_rows <= data_segment.n_cols - p_response + 1) {
    x_t_x = eye<mat>(
      data_segment.n_cols - p_response, data_segment.n_cols - p_response
    );
  } else {
    x_t_x = x.t() * x;
  }

  mat par = solve(x_t_x, x.t()) * y;
  mat residuals = y - x * par;
  double value =
    p_response * std::log(2.0 * M_PI) + log_det_sympd(variance_estimate);
  value *= data_segment.n_rows;
  value += trace(solve(variance_estimate, residuals.t() * residuals));
  return {{par}, {residuals}, value / 2};
}

CostResult Fastcpd::get_nll_variance(
  const unsigned int segment_start,
  const unsigned int segment_end
) {
  const unsigned int segment_length = segment_end - segment_start + 1;

  double det_value = det(arma::reshape(
    zero_data.row(segment_end + 1) - zero_data.row(segment_start), d, d
  ) / segment_length);
  if (segment_length < d) {
    unsigned int approximate_segment_start;
    unsigned int approximate_segment_end;
    if (segment_start >= d) {
      approximate_segment_start = segment_start - d;
    } else {
      approximate_segment_start = 0;
    }
    if (segment_end < data_n_rows - d) {
      approximate_segment_end = segment_end + d;
    } else {
      approximate_segment_end = data_n_rows - 1;
    }
    det_value = det(arma::reshape(
      zero_data.row(approximate_segment_end + 1) -
        zero_data.row(approximate_segment_start), d, d
    ) / (approximate_segment_end - approximate_segment_start + 1));
  }

  return {
    {zeros<mat>(d, d)},
    {mat()},
    (std::log(2.0 * M_PI) * d + d + log(det_value)) * segment_length / 2.0
  };
}

}  // namespace fastcpd::classes
