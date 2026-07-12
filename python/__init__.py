r"""
# What is fastcpd?

fastcpd is a Python package for fast change point detection in multivariate
time series data. It provides efficient algorithms to identify points in time
where the statistical properties of a sequence of observations change.

# Supported families

| Function | Family | Description |
|---|---|---|
| `segmentation.mean` | mean | Change in mean (univariate or multivariate) |
| `segmentation.variance` | variance | Change in variance (univariate or multivariate) |
| `segmentation.meanvariance` | meanvariance | Change in mean and/or variance |
| `segmentation.exponential` | exponential | Change in exponential rate parameter |
| `segmentation.var` | mgaussian | Change in VAR(p) model coefficients |
| `segmentation.lasso` | lasso | Change in penalised linear regression coefficients |
| `segmentation.garch` | garch | Change in GARCH(p,q) model parameters |
| `segmentation.lm` | gaussian | Change in linear regression coefficients |
| `segmentation.binomial` | binomial | Change in logistic regression coefficients |
| `segmentation.poisson` | poisson | Change in Poisson regression coefficients |
| `segmentation.arma` | arma | Change in ARMA(p,q) model parameters (q > 0) |
| `segmentation.ar` | ar | Change in AR(p) model parameters (OLS on lags) |
| `segmentation.arima` | arima | Change in ARIMA(p,d,q) model parameters |

Not supported: `custom` (requires Python/R callbacks). Pure MA models are
accessible via ``arima(data, order=(0, 0, q))``.

# Quickstart

``` python
import fastcpd.segmentation
from numpy import concatenate
from numpy.random import normal, multivariate_normal

# Mean change detection (multivariate)
covariance_mat = [[100, 0, 0], [0, 100, 0], [0, 0, 100]]
data = concatenate((multivariate_normal([0, 0, 0], covariance_mat, 300),
                    multivariate_normal([50, 50, 50], covariance_mat, 400),
                    multivariate_normal([2, 2, 2], covariance_mat, 300)))
result = fastcpd.segmentation.mean(data)
print(result.cp_set)  # [300, 700]

# Variance change detection (univariate)
data = concatenate((normal(0, 1, 500), normal(0, 5, 500)))
result = fastcpd.segmentation.variance(data)
print(result.cp_set)  # near 500

# Lasso regression change detection
from numpy import array
from numpy.random import randn, seed
import numpy as np
seed(7)
n, p = 400, 5
X = randn(n, p)
y = concatenate([X[:200] @ array([3., 0, 0, 0, 0]) + randn(200) * 0.1,
                 X[200:] @ array([0, 0, 0, 0, -3.]) + randn(200) * 0.1])
data = np.column_stack([y, X])
result = fastcpd.segmentation.lasso(data)
print(result.cp_set)  # [200]
```

"""

import os as _os
import sys as _sys

if _sys.platform == 'win32':
    # Python 3.8+ no longer searches PATH for DLL dependencies of extension
    # modules. Register the package directory so Windows finds libopenblas.dll,
    # which is bundled alongside interface.pyd by setup.py.
    _os.add_dll_directory(_os.path.dirname(_os.path.abspath(__file__)))

del _os, _sys

__version__ = "0.19.0"

from fastcpd.confidence import confint  # noqa: E402,F401
from fastcpd.segmentation import (  # noqa: E402
    CpdResult,
    ar,
    arima,
    arma,
    binomial,
    detect,
    detect_ar,
    detect_arima,
    detect_arma,
    detect_binomial,
    detect_exponential,
    detect_garch,
    detect_kcp,
    detect_kernel,
    detect_lasso,
    detect_linear_regression,
    detect_lm,
    detect_logistic_regression,
    detect_mean,
    detect_mean_variance,
    detect_meanvariance,
    detect_poisson,
    detect_poisson_regression,
    detect_quantile,
    detect_quantile_regression,
    detect_rank,
    detect_time_series,
    detect_ts,
    detect_var,
    detect_variance,
    exponential,
    garch,
    kernel,
    kcp,
    lasso,
    lm,
    mean,
    meanvariance,
    poisson,
    quantile,
    rank,
    time_series,
    ts,
    var,
    variance,
)
from fastcpd.variance_estimation import (  # noqa: E402
    VarianceArmaResult,
    estimate_variance,
    estimate_variance_arma,
    estimate_variance_linear_regression,
    estimate_variance_lm,
    estimate_variance_mean,
    estimate_variance_median,
)

__all__ = [
    'CpdResult',
    'VarianceArmaResult',
    'ar',
    'arima',
    'arma',
    'binomial',
    'confint',
    'detect',
    'detect_ar',
    'detect_arima',
    'detect_arma',
    'detect_binomial',
    'detect_exponential',
    'detect_garch',
    'detect_kcp',
    'detect_kernel',
    'detect_lasso',
    'detect_linear_regression',
    'detect_lm',
    'detect_logistic_regression',
    'detect_mean',
    'detect_mean_variance',
    'detect_meanvariance',
    'detect_poisson',
    'detect_poisson_regression',
    'detect_quantile',
    'detect_quantile_regression',
    'detect_rank',
    'detect_time_series',
    'detect_ts',
    'detect_var',
    'detect_variance',
    'exponential',
    'estimate_variance',
    'estimate_variance_arma',
    'estimate_variance_linear_regression',
    'estimate_variance_lm',
    'estimate_variance_mean',
    'estimate_variance_median',
    'garch',
    'kernel',
    'kcp',
    'lasso',
    'lm',
    'mean',
    'meanvariance',
    'poisson',
    'quantile',
    'rank',
    'time_series',
    'ts',
    'var',
    'variance',
]
