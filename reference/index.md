# Package index

## Main function

All implementation of fastcpd is unified into one single function
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md).

- [`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) : Find
  change points efficiently

## Wrapper functions

Simplified functions are provided for each family.

### Time series

[`fastcpd_ts()`](https://fastcpd.xingchi.li/reference/fastcpd_ts.md) and
[`fastcpd.ts()`](https://fastcpd.xingchi.li/reference/fastcpd_ts.md) are
the main functions for time series data (also wrapper functions of
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md)).

- [`fastcpd_ar()`](https://fastcpd.xingchi.li/reference/fastcpd_ar.md)
  [`fastcpd.ar()`](https://fastcpd.xingchi.li/reference/fastcpd_ar.md) :
  Find change points efficiently in AR(\\p\\) models
- [`fastcpd_arima()`](https://fastcpd.xingchi.li/reference/fastcpd_arima.md)
  [`fastcpd.arima()`](https://fastcpd.xingchi.li/reference/fastcpd_arima.md)
  : Find change points efficiently in ARIMA(\\p\\, \\d\\, \\q\\) models
- [`fastcpd_arma()`](https://fastcpd.xingchi.li/reference/fastcpd_arma.md)
  [`fastcpd.arma()`](https://fastcpd.xingchi.li/reference/fastcpd_arma.md)
  : Find change points efficiently in ARMA(\\p\\, \\q\\) models
- [`fastcpd_garch()`](https://fastcpd.xingchi.li/reference/fastcpd_garch.md)
  [`fastcpd.garch()`](https://fastcpd.xingchi.li/reference/fastcpd_garch.md)
  : Find change points efficiently in GARCH(\\p\\, \\q\\) models
- [`fastcpd_ts()`](https://fastcpd.xingchi.li/reference/fastcpd_ts.md)
  [`fastcpd.ts()`](https://fastcpd.xingchi.li/reference/fastcpd_ts.md) :
  Find change points efficiently in time series data
- [`fastcpd_var()`](https://fastcpd.xingchi.li/reference/fastcpd_var.md)
  [`fastcpd.var()`](https://fastcpd.xingchi.li/reference/fastcpd_var.md)
  : Find change points efficiently in VAR(\\p\\) models

### Unlabeled data

Used for data without response variables, for example, mean change and
variance change.

- [`fastcpd_mean()`](https://fastcpd.xingchi.li/reference/fastcpd_mean.md)
  [`fastcpd.mean()`](https://fastcpd.xingchi.li/reference/fastcpd_mean.md)
  : Find change points efficiently in mean change models
- [`fastcpd_variance()`](https://fastcpd.xingchi.li/reference/fastcpd_variance.md)
  [`fastcpd.variance()`](https://fastcpd.xingchi.li/reference/fastcpd_variance.md)
  : Find change points efficiently in variance change models
- [`fastcpd_meanvariance()`](https://fastcpd.xingchi.li/reference/fastcpd_meanvariance.md)
  [`fastcpd.meanvariance()`](https://fastcpd.xingchi.li/reference/fastcpd_meanvariance.md)
  [`fastcpd_mv()`](https://fastcpd.xingchi.li/reference/fastcpd_meanvariance.md)
  [`fastcpd.mv()`](https://fastcpd.xingchi.li/reference/fastcpd_meanvariance.md)
  : Find change points efficiently in mean variance change models

### Regression data

Detect change points in the coefficients of regression-type data.

- [`fastcpd_binomial()`](https://fastcpd.xingchi.li/reference/fastcpd_binomial.md)
  [`fastcpd.binomial()`](https://fastcpd.xingchi.li/reference/fastcpd_binomial.md)
  : Find change points efficiently in logistic regression models
- [`fastcpd_lasso()`](https://fastcpd.xingchi.li/reference/fastcpd_lasso.md)
  [`fastcpd.lasso()`](https://fastcpd.xingchi.li/reference/fastcpd_lasso.md)
  : Find change points efficiently in penalized linear regression models
- [`fastcpd_lm()`](https://fastcpd.xingchi.li/reference/fastcpd_lm.md)
  [`fastcpd.lm()`](https://fastcpd.xingchi.li/reference/fastcpd_lm.md) :
  Find change points efficiently in linear regression models
- [`fastcpd_poisson()`](https://fastcpd.xingchi.li/reference/fastcpd_poisson.md)
  [`fastcpd.poisson()`](https://fastcpd.xingchi.li/reference/fastcpd_poisson.md)
  : Find change points efficiently in Poisson regression models

## Utility functions

The following functions help with visualization and analyzation of the
data.

### Variance estimation

- [`variance_arma()`](https://fastcpd.xingchi.li/reference/variance_arma.md)
  [`variance.arma()`](https://fastcpd.xingchi.li/reference/variance_arma.md)
  : Variance estimation for ARMA model with change points
- [`variance_lm()`](https://fastcpd.xingchi.li/reference/variance_lm.md)
  [`variance.lm()`](https://fastcpd.xingchi.li/reference/variance_lm.md)
  : Variance estimation for linear models with change points
- [`variance_mean()`](https://fastcpd.xingchi.li/reference/variance_mean.md)
  [`variance.mean()`](https://fastcpd.xingchi.li/reference/variance_mean.md)
  : Variance estimation for mean change models
- [`variance_median()`](https://fastcpd.xingchi.li/reference/variance_median.md)
  [`variance.median()`](https://fastcpd.xingchi.li/reference/variance_median.md)
  : Variance estimation for median change models

### Class methods

- [`plot(`*`<fastcpd>`*`)`](https://fastcpd.xingchi.li/reference/plot.md)
  [`plot(`*`<fastcpd>`*`,`*`<missing>`*`)`](https://fastcpd.xingchi.li/reference/plot.md)
  : Plot the data and the change points for a fastcpd object
- [`print(`*`<fastcpd>`*`)`](https://fastcpd.xingchi.li/reference/print.md)
  : Print the call and the change points for a fastcpd object
- [`show(`*`<fastcpd>`*`)`](https://fastcpd.xingchi.li/reference/show.md)
  : Show the available methods for a fastcpd object
- [`summary(`*`<fastcpd>`*`)`](https://fastcpd.xingchi.li/reference/summary.md)
  : Show the summary of a fastcpd object

## Data

fastcpd comes with a selection of built-in datasets that are used in
examples to illustrate various change point detection challenges.

- [`bitcoin`](https://fastcpd.xingchi.li/reference/bitcoin.md) : Bitcoin
  Market Price (USD)
- [`occupancy`](https://fastcpd.xingchi.li/reference/occupancy.md) :
  Occupancy Detection Data Set
- [`transcriptome`](https://fastcpd.xingchi.li/reference/transcriptome.md)
  : Transcription Profiling of 57 Human Bladder Carcinoma Samples
- [`uk_seatbelts`](https://fastcpd.xingchi.li/reference/uk_seatbelts.md)
  : UK Seatbelts Data
- [`well_log`](https://fastcpd.xingchi.li/reference/well_log.md) :
  Well-log Dataset from Numerical Bayesian Methods Applied to Signal
  Processing

## Main class

The results of
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) will be
stored in the class `fastcpd` with accompanied several utility
functions.

- [`fastcpd-class`](https://fastcpd.xingchi.li/reference/fastcpd-class.md)
  :

  An S4 class to store the output created with
  [`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md)
