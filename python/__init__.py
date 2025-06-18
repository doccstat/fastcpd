r"""
# What is fastcpd?

fastcpd is a Python package for fast change point detection in multivariate
time series data. It provides efficient algorithms to identify points in time
where the statistical properties of a sequence of observations change.

# Quickstart

``` python
import fastcpd.segmentation
from numpy import concatenate
from numpy.random import normal, multivariate_normal
covariance_mat = [[100, 0, 0], [0, 100, 0], [0, 0, 100]]
data = concatenate((multivariate_normal([0, 0, 0], covariance_mat, 300),
                    multivariate_normal([50, 50, 50], covariance_mat, 400),
                    multivariate_normal([2, 2, 2], covariance_mat, 300)))
fastcpd.segmentation.mean(data)
```

"""

__version__ = "0.18.0"
