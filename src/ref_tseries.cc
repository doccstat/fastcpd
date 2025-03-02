// This is a modified copy of tseries/R/garch.R

#include "ref_tseries.h"

#include <algorithm>
#include <cmath>
#include <limits>

using ::Rcpp::CharacterVector;
using ::Rcpp::Function;
using ::Rcpp::NumericMatrix;

using ::Rcpp::as;
using ::Rcpp::Named;
using ::Rcpp::sqrt;
using ::Rcpp::stop;
using ::Rcpp::warning;
using ::Rcpp::wrap;
using ::std::max;
using ::std::pow;
using ::std::to_string;

// Declare the C functions with C linkage.
extern "C" {
void tseries_fit_garch(double *y, int *n, double *par, int *p, int *q,
                       int *itmax, double *afctol, double *rfctol,
                       double *xctol, double *xftol, double *fret, int *agrad,
                       int *trace);

void tseries_pred_garch(double *y, double *h, int *n, double *par, int *p,
                        int *q, int *genuine);

void tseries_ophess_garch(double *y, int *n, double *par, double *he, int *p,
                          int *q);
}

List garch(const colvec &x_, const colvec &order_, Nullable<string> series,
           int maxiter, bool trace, Nullable<NumericVector> start, string grad,
           Nullable<double> abstol_, Nullable<double> reltol_,
           Nullable<double> xtol_, Nullable<double> falsetol_) {
  NumericVector x = wrap(x_);
  IntegerVector order = wrap(order_);
  // Set machine-epsilon based defaults:
  double eps = std::numeric_limits<double>::epsilon();
  double abstol =
      abstol_.isNull() ? max(1e-20, pow(eps, 2.0)) : as<double>(abstol_);
  double reltol =
      reltol_.isNull() ? max(1e-10, pow(eps, 2.0 / 3.0)) : as<double>(reltol_);
  double xtol = xtol_.isNull() ? sqrt(eps) : as<double>(xtol_);
  double falsetol = falsetol_.isNull() ? 1e2 * eps : as<double>(falsetol_);

  // Set up gradient evaluation method.
  int agrad;
  if (grad == "analytical") {
    agrad = 1;
  } else if (grad == "numerical") {
    agrad = 0;
  } else {
    stop("grad must be either 'analytical' or 'numerical'");
  }

  // Determine the series name.
  string series_name = series.isNotNull() ? as<string>(series) : "x";

  // Check for NA’s in x.
  for (int i = 0; i < x.size(); i++) {
    if (NumericVector::is_na(x[i])) stop("NAs in x");
  }

  int n = x.size();

  // 'order' must be of length 2.
  if (order.size() != 2) stop("order must be of length 2");
  // In our convention, order[0] is p and order[1] is q.
  int p = order[0];
  int q = order[1];
  int ncoef = p + q + 1;

  double small = 0.05;

  // Set up starting coefficients.
  NumericVector coef;
  if (start.isNotNull()) {
    coef = as<NumericVector>(start);
  } else {
    // Default: first coefficient a0 is var(x)*(1 - small*(ncoef-1))
    double sum = 0.0, sumsq = 0.0;
    for (int i = 0; i < n; i++) {
      sum += x[i];
      sumsq += x[i] * x[i];
    }
    double mean = sum / n;
    double var = (sumsq - n * mean * mean) / (n - 1);
    double a0 = var * (1.0 - small * (ncoef - 1));
    coef = NumericVector(ncoef, small);
    coef[0] = a0;
  }
  if (coef.size() != ncoef) stop("incorrect length of coef");

  // Prepare trace flag as integer.
  int itrace = trace ? 1 : 0;

  // --- Call the C functions directly ---
  // 1. Fit the model.
  double nlikeli = 1e+10;
  tseries_fit_garch(x.begin(), &n, coef.begin(), &p, &q, &maxiter, &abstol,
                    &reltol, &xtol, &falsetol, &nlikeli, &agrad, &itrace);

  // 2. Compute predictions.
  NumericVector pred(n);
  int genuine_int = 0;  // FALSE
  tseries_pred_garch(x.begin(), pred.begin(), &n, coef.begin(), &p, &q,
                     &genuine_int);

  // 3. Compute the observed Hessian.
  NumericMatrix com_hess(ncoef, ncoef);
  tseries_ophess_garch(x.begin(), &n, coef.begin(), com_hess.begin(), &p, &q);

  // --- Compute the inverse (vcov) using QR decomposition via R's functions ---
  Function qr("qr");
  List qr_list = qr(com_hess);
  int rank = as<int>(qr_list["rank"]);

  NumericMatrix vc(ncoef, ncoef);
  if (rank != ncoef) {
    // If the Hessian is singular, fill vc with NA’s and warn.
    for (int i = 0; i < ncoef; i++)
      for (int j = 0; j < ncoef; j++) vc(i, j) = NA_REAL;
    warning("singular information");
  } else {
    Function solve("solve");
    vc = solve(com_hess);
  }

  // --- Post-processing: create fitted values and residuals ---
  // sigt = sqrt(pred)
  NumericVector sigt(n);
  for (int i = 0; i < n; i++) sigt[i] = sqrt(pred[i]);

  // For the first max(p, q) observations, set sigt to NA.
  int max_order = max(p, q);
  for (int i = 0; i < max_order && i < n; i++) sigt[i] = NA_REAL;

  // Create fitted values matrix with two columns: sigt and -sigt.
  NumericMatrix f(n, 2);
  for (int i = 0; i < n; i++) {
    f(i, 0) = sigt[i];
    f(i, 1) = -sigt[i];
  }
  CharacterVector fcolnames = CharacterVector::create("sigt", "-sigt");
  f.attr("dimnames") = List::create(R_NilValue, fcolnames);

  // Compute standardized residuals e = x / sigt.
  NumericVector e(n);
  for (int i = 0; i < n; i++) {
    if (NumericVector::is_na(sigt[i]))
      e[i] = NA_REAL;
    else
      e[i] = x[i] / sigt[i];
  }

  // If x is a time series (has a "tsp" attribute), propagate it.
  if (x.hasAttribute("tsp")) {
    SEXP tsp_attr = x.attr("tsp");
    f.attr("tsp") = tsp_attr;
    e.attr("tsp") = tsp_attr;
    f.attr("class") = "ts";
    e.attr("class") = "ts";
  }

  // --- Naming the results ---
  // Name the order vector.
  CharacterVector order_names = CharacterVector::create("p", "q");
  order.attr("names") = order_names;

  // Create coefficient names: first "a0", then "a1",..., for q, then "b1",...
  // for p.
  CharacterVector coef_names(ncoef);
  coef_names[0] = "a0";
  int idx = 1;
  for (int i = 0; i < q; i++) coef_names[idx++] = "a" + to_string(i + 1);
  for (int i = 0; i < p; i++) coef_names[idx++] = "b" + to_string(i + 1);
  coef.attr("names") = coef_names;
  vc.attr("dimnames") = List::create(coef_names, coef_names);

  // Get frequency from the time series attributes (if present).
  double xfreq = NA_REAL;
  if (x.hasAttribute("tsp")) {
    NumericVector tsp_attr = x.attr("tsp");
    if (tsp_attr.size() >= 3) xfreq = tsp_attr[2];
  }

  // --- Assemble the output list ---
  List result = List::create(
      Named("order") = order, Named("coef") = coef, Named("n.likeli") = nlikeli,
      Named("n.used") = n, Named("residuals") = e, Named("fitted.values") = f,
      Named("series") = series_name, Named("frequency") = xfreq,
      Named("call") = R_NilValue, Named("vcov") = vc);
  result.attr("class") = "garch";

  return result;
}
