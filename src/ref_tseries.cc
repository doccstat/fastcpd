// garch_wrappers.cpp
#include <Rcpp.h>

// Declare the C functions with C linkage.
extern "C"
{
  void tseries_fit_garch(double *y,
                         int *n,
                         double *par,
                         int *p,
                         int *q,
                         int *itmax,
                         double *afctol,
                         double *rfctol,
                         double *xctol,
                         double *xftol,
                         double *fret,
                         int *agrad,
                         int *trace);

  void tseries_pred_garch(double *y,
                          double *h,
                          int *n,
                          double *par,
                          int *p,
                          int *q,
                          int *genuine);

  void tseries_ophess_garch(double *y,
                            int *n,
                            double *par,
                            double *he,
                            int *p,
                            int *q);
}

// [[Rcpp::export]]
Rcpp::List fit_garch_wrapper(Rcpp::NumericVector x,
                             int n,
                             Rcpp::NumericVector coef,
                             int p,
                             int q,
                             int maxiter,
                             double abstol,
                             double reltol,
                             double xtol,
                             double falsetol,
                             int agrad,
                             int trace)
{
  // Initialize likelihood output.
  double nlikeli = 1e+10;

  // Call the C function.
  tseries_fit_garch(x.begin(), &n, coef.begin(),
                    &p, &q,
                    &maxiter,
                    &abstol, &reltol, &xtol, &falsetol,
                    &nlikeli, &agrad, &trace);

  // Return the updated coefficients and likelihood as a list.
  return Rcpp::List::create(Rcpp::Named("coef") = coef,
                            Rcpp::Named("nlikeli") = nlikeli);
}

// [[Rcpp::export]]
Rcpp::NumericVector pred_garch_wrapper(Rcpp::NumericVector x,
                                       int n,
                                       Rcpp::NumericVector coef,
                                       int p,
                                       int q,
                                       bool genuine)
{
  // Allocate a numeric vector for the predictions.
  Rcpp::NumericVector e(n);

  // Convert the boolean 'genuine' to an integer (1 for true, 0 for false)
  int genuine_int = genuine ? 1 : 0;

  // Call the C function.
  tseries_pred_garch(x.begin(), e.begin(), &n, coef.begin(),
                     &p, &q, &genuine_int);
  return e;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix ophess_garch_wrapper(Rcpp::NumericVector x,
                                         int n,
                                         Rcpp::NumericVector coef,
                                         int p,
                                         int q)
{
  // Typically, the number of parameters in a GARCH model is p+q+1.
  int npar = p + q + 1;

  // Preallocate a matrix for the Hessian.
  Rcpp::NumericMatrix hess(npar, npar);

  // Call the C function.
  tseries_ophess_garch(x.begin(), &n, coef.begin(), hess.begin(),
                       &p, &q);
  return hess;
}
