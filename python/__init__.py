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

Families not yet supported in the Python binding (require R packages at
runtime): `arma`, `arima`, `ma`, `garch`, `lm`, `binomial`, `poisson`,
`custom`.

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

__version__ = "0.19.0"
