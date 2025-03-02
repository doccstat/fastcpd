// This is a modified copy of tseries/src/garch.c

/* Copyright (C) 1997-1999  Adrian Trapletti

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public
   License as published by the Free Software Foundation; either
   version 2 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public
   License along with this library; if not, write to the Free
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

   GARCH estimation

   Reference: T. Bollerslev (1986): Generalized Autoregressive Conditional
   Heteroscedasticity, Journal of Econometrics 31, 307-327. */

#include <R.h>

extern void F77_NAME(dsumsl)(int *n, double *d, double *x,
                             void (*)(int *, double *, int *, double *, int *,
                                      double *, void (*)(void)),
                             void (*)(int *, double *, int *, double *, int *,
                                      double *, void (*)(void)),
                             int *iv, int *liv, int *lv, double *v, int *uiparm,
                             double *urparm, void (*)(void));
extern void F77_NAME(dsmsno)(int *n, double *d, double *x,
                             void (*)(int *, double *, int *, double *, int *,
                                      double *, void (*)(void)),
                             int *iv, int *liv, int *lv, double *v, int *uiparm,
                             double *urparm, void (*)(void));
extern void F77_NAME(ddeflt)(int *alg, int *iv, int *liv, int *lv, double *v);

#define BIG 1.0e+10 /* function value if the parameters are invalid */

static double dsqrarg;
#define DSQR(a) ((dsqrarg = (a)) == 0.0 ? 0.0 : dsqrarg * dsqrarg)

static double dmaxarg1, dmaxarg2;
#define DMAX(a, b)                 \
  (dmaxarg1 = (a), dmaxarg2 = (b), \
   (dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))

struct garch_handler /* used to set up the additional parameters used in calcf
                        and calcg */
{
  double *y;  /* the time series to fit */
  double *h;  /* the conditional variance (cv) */
  double *dh; /* dh_i/dp_j */
  int n;      /* the number of observations */
  int p, q;   /* GARCH(p,q) */
};

static struct garch_handler garch_h;

static void F77_SUB(calcf)(int *pq, double *p, int *nf, double *f, int *uiparm,
                           double *urparm, void (*F77_SUB(ufparm))(void))
/* compute negative log likelihood apart from the constant and the pre-sample
   values */
{
  int i, j, ok;
  int maxpq = (int)DMAX(garch_h.p, garch_h.q);
  double temp = 0.0;
  double sum = 0.0;

  ok = 1;
  if (p[0] <= 0.0) ok = 0;
  for (i = 1; i < (*pq); i++)
    if (p[i] < 0.0) ok = 0;
  if (ok) /* parameters are valid */
  {
    for (i = maxpq; i < garch_h.n; i++) /* loop over time */
    {                                   /* compute cv at time i */
      temp = p[0];                      /* compute ARCH part of cv */
      for (j = 1; j <= garch_h.q; j++) temp += p[j] * DSQR(garch_h.y[i - j]);
      for (j = 1; j <= garch_h.p; j++) /* compute GARCH part of cv */
        temp += p[garch_h.q + j] * garch_h.h[i - j];
      sum += log(temp) + DSQR(garch_h.y[i]) / temp; /* compute eq. 18 */
      garch_h.h[i] = temp;                          /* assign cv at time i */
    }
    (*f) = 0.5 * sum;
  } else /* parameters are invalid */
    (*f) = BIG;
}

static void F77_SUB(calcg)(int *pq, double *p, int *nf, double *dp, int *uiparm,
                           double *urparm, void (*F77_SUB(ufparm))(void))
/* compute derivative of negative log likelihood */
{
  int i, j, k;
  int maxpq = (int)DMAX(garch_h.p, garch_h.q);
  double temp1, temp2, temp3;

  for (k = 0; k < (*pq); k++) /* initialize */
    dp[k] = 0.0;
  for (i = maxpq; i < garch_h.n; i++) /* loop over time */
  {               /* compute cv at time i and derivatives dh_i/dp_j */
    temp1 = p[0]; /* compute ARCH part of cv */
    for (j = 1; j <= garch_h.q; j++) temp1 += p[j] * DSQR(garch_h.y[i - j]);
    for (j = 1; j <= garch_h.p; j++) /* compute GARCH part of cv */
      temp1 += p[garch_h.q + j] * garch_h.h[i - j];
    garch_h.h[i] = temp1; /* assign cv at time i */
    temp2 = 0.5 * (1.0 - DSQR(garch_h.y[i]) / temp1) /
            temp1; /* compute dl_i/dh_i, eq. 19 */
    temp3 = 1.0;   /* compute dh_i/dp_0, eq. 21 */
    for (j = 1; j <= garch_h.p; j++)
      temp3 += p[garch_h.q + j] * garch_h.dh[(*pq) * (i - j)];
    garch_h.dh[(*pq) * i] = temp3; /* assign dh_i/dp_0 */
    dp[0] += temp2 * temp3; /* assign dl_i/dp_0 = dl_i/dh_i * dh_i/dp_0 */
    for (k = 1; k <= garch_h.q;
         k++) /* compute dl_i/dp_k for the ARCH part, eq. 19 */
    {
      temp3 = DSQR(garch_h.y[i - k]); /* compute dh_i/dp_k, eq. 21 */
      for (j = 1; j <= garch_h.p; j++)
        temp3 += p[garch_h.q + j] * garch_h.dh[(*pq) * (i - j) + k];
      garch_h.dh[(*pq) * i + k] = temp3; /* assign dh_i/dp_k */
      dp[k] += temp2 * temp3; /* assign dl_i/dp_k = dl_i/dh_i * dh_i/dp_k */
    }
    for (k = 1; k <= garch_h.p;
         k++) /* compute dl_i/dp_k for the GARCH part, eq. 19 */
    {
      temp3 = garch_h.h[i - k]; /* compute dh_i/dp_k, eq. 21 */
      for (j = 1; j <= garch_h.p; j++)
        temp3 += p[garch_h.q + j] * garch_h.dh[(*pq) * (i - j) + garch_h.q + k];
      garch_h.dh[(*pq) * i + garch_h.q + k] = temp3; /* assign dh_i/dp_k */
      dp[garch_h.q + k] +=
          temp2 * temp3; /* assign dl_i/dp_k = dl_i/dh_i * dh_i/dp_k */
    }
  }
}

static void F77_SUB(ufparm)(void) { error("fatal error in fit_garch ()\n"); }

void tseries_fit_garch(double *y, int *n, double *par, int *p, int *q,
                       int *itmax, double *afctol, double *rfctol,
                       double *xctol, double *xftol, double *fret, int *agrad,
                       int *trace)
/* fit a GARCH (p, q) model

Input:

y[0..n-1]      series to fit
par[0..p+q]    initial parameter estimates
p, q           model orders
itmax          maximum number of iterations
afctol         absolute function convergence tolerance
   rfctol         relative function convergence tolerance
   xctol          x-convergence tolerance
   xftol          false convergence tolerance
agrad          estimation with analytical/numerical gradient
trace          output yes/no

Output:

par[0..p+q]    parameter estimates at minimum
fret           function value at minimum
*/
{
  int i, j, pq, liv, lv, alg;
  double *d, *v;
  int *iv;
  double var;

  /* set up general optimizer parameters to default values */
  pq = (*p) + (*q) + 1;
  d = R_Calloc(pq, double);
  for (i = 0; i < pq; i++) d[i] = 1.0;
  liv = 60;
  iv = R_Calloc(liv, int);
  lv = 77 + pq * (pq + 17) / 2;
  v = R_Calloc(lv, double);
  alg = 2;
  F77_CALL(ddeflt)(&alg, iv, &liv, &lv, v);
  iv[0] = 12;

  /* set up user defined optimizer parameters */
  iv[16] = 2 * (*itmax);
  iv[17] = (*itmax);
  if (*trace)
    iv[20] = 6;
  else
    iv[20] = 0;
  v[30] = (*afctol);
  v[31] = (*rfctol);
  v[32] = (*xctol);
  v[33] = (*xftol);

  /* set handler values */
  garch_h.p = (*p);
  garch_h.q = (*q);
  garch_h.n = (*n);
  garch_h.y = y;
  garch_h.h = R_Calloc((*n), double);
  garch_h.dh = R_Calloc((*n) * pq, double);
  var = 0.0;
  for (i = 0; i < (*n); i++) /* estimate unconditional variance (uv) */
    var += DSQR(y[i]);
  var /= (double)(*n);
  for (i = 0; i < (int)DMAX((*p), (*q)); i++) /* initialize */
  {
    garch_h.h[i] = var;       /* with uv */
    garch_h.dh[pq * i] = 1.0; /* dh_i/dp_0 with 1 */
    for (j = 1; j < pq; j++)
      garch_h.dh[pq * i + j] = 0.0; /* dh_i/dp_j with 0 */
  }

  if (*agrad) /* estimation with analytical gradient */
  {
    if (*trace)
      Rprintf("\n ***** ESTIMATION WITH ANALYTICAL GRADIENT ***** \n\n");
    F77_CALL(dsumsl)(&pq, d, par, F77_SUB(calcf), F77_SUB(calcg), iv, &liv, &lv,
                     v, NULL, NULL, F77_SUB(ufparm));
    if (*trace) Rprintf("\n");
  } else /* estimation with numerical gradient */
  {
    if (*trace)
      Rprintf("\n ***** ESTIMATION WITH NUMERICAL GRADIENT ***** \n\n");
    F77_CALL(dsmsno)(&pq, d, par, F77_SUB(calcf), iv, &liv, &lv, v, NULL, NULL,
                     F77_SUB(ufparm));
    if (*trace) Rprintf("\n");
  }

  /* return function value */
  (*fret) = v[9];

  /* free memory */
  R_Free(d);
  R_Free(iv);
  R_Free(v);
  R_Free(garch_h.h);
  R_Free(garch_h.dh);
}

void tseries_pred_garch(double *y, double *h, int *n, double *par, int *p,
                        int *q, int *genuine)
/* predict cv with a GARCH (p, q) model

Input:

y[0..n-1]    series to predict
par[0..p+q]  parameters of the GARCH (p, q)
p, q         model orders
genuine      logical indicating if a genuine prediction is computed

Output:

h[0..N]      predicted cv, where N = n for genuine prediction, and N = n-1
otherwise
*/
{
  double var, temp;
  int i, j, maxpq, N;

  if (*genuine)
    N = (*n) + 1;
  else
    N = (*n);
  maxpq = (int)DMAX((*p), (*q));
  var = 0.0;
  for (i = 1; i <= (*p) + (*q); i++) /* compute uv */
    var += par[i];
  var = par[0] / (1.0 - var);
  for (i = 0; i < maxpq; i++) /* initialize with uv */
    h[i] = var;
  for (i = maxpq; i < N; i++) /* loop over time */
  {                           /* compute cv at time i */
    temp = par[0];            /* compute ARCH part of cv */
    for (j = 1; j <= (*q); j++) temp += par[j] * DSQR(y[i - j]);
    for (j = 1; j <= (*p); j++) /* compute GARCH part of cv */
      temp += par[(*q) + j] * h[i - j];
    h[i] = temp; /* assign cv at time i */
  }
}

void tseries_ophess_garch(double *y, int *n, double *par, double *he, int *p,
                          int *q)
/* Compute outer product approximation of the hessian of the
negative log likelihood of a GARCH (p, q) model at given parameter
estimates

Input:

y[0..n-1]              time series
par[0..p+q]            parameter estimates of the GARCH (p, q)
p, q                   model orders

Output:

he[0..(p+q+1)*(p+q+1)-1]      predicted cv
*/
{
  double var, temp1, temp2, temp3;
  int i, j, k, pq;
  double *h, *dh, *dpar;

  pq = (*p) + (*q) + 1;
  h = R_Calloc((*n), double);
  dh = R_Calloc((*n) * pq, double);
  dpar = R_Calloc(pq, double);
  var = 0.0;
  for (i = 0; i < (*n); i++) /* estimate uv */
    var += DSQR(y[i]);
  var /= (double)(*n);
  for (i = 0; i < (int)DMAX((*p), (*q)); i++) /* initialize */
  {
    h[i] = var;                                    /* with uv */
    dh[pq * i] = 1.0;                              /* dh_i/dp_0 with 1 */
    for (j = 1; j < pq; j++) dh[pq * i + j] = 0.0; /* dh_i/dp_j with 0 */
  }
  for (k = 0; k < pq; k++) /* initialize */
    for (j = 0; j < pq; j++) he[pq * k + j] = 0.0;
  for (i = (int)DMAX((*p), (*q)); i < (*n); i++) /* loop over time */
  {                 /* compute cv at time i and derivatives dh_i/dp_j */
    temp1 = par[0]; /* compute ARCH part of cv */
    for (j = 1; j <= (*q); j++) temp1 += par[j] * DSQR(y[i - j]);
    for (j = 1; j <= (*p); j++) /* compute GARCH part of cv */
      temp1 += par[(*q) + j] * h[i - j];
    h[i] = temp1; /* assign cv at time i */
    temp2 = 0.5 * (1.0 - DSQR(y[i]) / temp1) /
            temp1; /* compute dl_i/dh_i, eq. 19 */
    temp3 = 1.0;   /* compute dh_i/dp_0, eq. 21 */
    for (j = 1; j <= (*p); j++) temp3 += par[(*q) + j] * dh[pq * (i - j)];
    dh[pq * i] = temp3;      /* assign dh_i/dp_0 */
    dpar[0] = temp2 * temp3; /* assign dl_i/dp_0 = dl_i/dh_i * dh_i/dp_0 */
    for (k = 1; k <= (*q);
         k++) /* compute dl_i/dp_k for the ARCH part, eq. 19 */
    {
      temp3 = DSQR(y[i - k]); /* compute dh_i/dp_k, eq. 21 */
      for (j = 1; j <= (*p); j++) temp3 += par[(*q) + j] * dh[pq * (i - j) + k];
      dh[pq * i + k] = temp3;  /* assign dh_i/dp_k */
      dpar[k] = temp2 * temp3; /* assign dl_i/dp_k = dl_i/dh_i * dh_i/dp_k */
    }
    for (k = 1; k <= (*p);
         k++) /* compute dl_i/dp_k for the GARCH part, eq. 19 */
    {
      temp3 = h[i - k]; /* compute dh_i/dp_k, eq. 21 */
      for (j = 1; j <= (*p); j++)
        temp3 += par[(*q) + j] * dh[pq * (i - j) + (*q) + k];
      dh[pq * i + (*q) + k] = temp3; /* assign dh_i/dp_k */
      dpar[(*q) + k] =
          temp2 * temp3; /* assign dl_i/dp_k = dl_i/dh_i * dh_i/dp_k */
    }
    for (k = 0; k < pq; k++) /* compute outer product approximation, p. 317 */
      for (j = 0; j < pq; j++) he[pq * k + j] += dpar[k] * dpar[j];
  }
  R_Free(h);
  R_Free(dh);
  R_Free(dpar);
}
