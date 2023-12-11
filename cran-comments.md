## Updates since last CRAN release (0.9.0)

### fastcpd 0.9.8

*   Preliminary support for ARMA(p, q) model with parameter `order = c(p, q)`
    and family `"arma"`.
*   Add `fastcpd.arma` / `fastcpd_arma` for ARMA(p, q) model.
*   Add adaptive increasing `beta` values.

### fastcpd 0.9.7

*   Add `lower` and `upper` parameters to denote the lower and upper bounds of
    the parameters.
*   Add line search.
*   Add hardcoded ARMA(3, 2).

### fastcpd 0.9.6

*   Add `bitcoin` and `well_log` data.
*   Add residual calculation for mean family.
*   Add plots for bitcoin data.
*   Fix residual calculation for time series data when the seegments are
    too small.
*   Handle variance estimation errors.

### fastcpd 0.9.5

*   Add wrapper functions of
    AR(p) family: `fastcpd.ar` / `fastcpd_ar`,
    ARIMA(p, d, q) family: `fastcpd.arima` / `fastcpd_arima`,
    GARCH(p, q) family: `fastcpd.garch` / `fastcpd_garch`,
    linear regression family: `fastcpd.lm` / `fastcpd_lm`,
    logistic regression family: `fastcpd.binomial` / `fastcpd_binomial`,
    poisson regression family: `fastcpd.poisson` / `fastcpd_poisson`,
    penalized linear regression family: `fastcpd.lasso` / `fastcpd_lasso`,
    MA(q) model: `fastcpd.ma` / `fastcpd_ma`,
    mean change: `fastcpd.mean` / `fastcpd_mean`,
    variance change: `fastcpd.variance` / `fastcpd_variance`,
    mean or variance change: `fastcpd.meanvariance` / `fastcpd_meanvariance` /
    `fastcpd.mv` / `fastcpd_mv`.
*   Replace `"gaussian"` family with `"lm"`.
*   Add progress bar.
*   Fix design matrix from formula bug.

### fastcpd 0.9.4

*   Fix sanity check.
*   Add small AR(1) data in gallery.
*   Fix VAR(p) model bug.
*   Add VAR(2) example in Gallery.
*   Remove commented code.

### fastcpd 0.9.3

*   Deprecate "vanilla" family by `vanilla_percentage` parameter.
*   Add check utility functions.
*   Add MA(4) example.
*   Fix the bug when `beta` is updated but the old `beta` is still in use.
*   Remove tests estimating the variance in the "gaussian" family dynamically.
*   Merge `beta` updating into `get_segment_statistics`.

### fastcpd 0.9.2

*   Add gallery to vignettes.
*   Remove cheatsheets pdf from the package.
*   Use `forecast` package for ARIMA model.
*   Use `fGarch` package for GARCH model.

### fastcpd 0.9.1

*   Wrap `&&` around `||` by parentheses.
*   Add ma(4) example using custom cost function.
*   Add full support for AR(p), MA(q) and ARIMA(p, d, q) models.
