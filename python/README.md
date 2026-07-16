# fastcpd for Python

`fastcpd` provides fast change-point detection for Python through the same
canonical C++ implementation used by the standalone `fastcpd-cpp` package.
The Python package is independently buildable: it does not require R, invoke
R code generation, or download another fastcpd repository.

## Install

```shell
python -m pip install fastcpd
```

```python
import numpy as np
from fastcpd import detect_mean

data = np.concatenate([np.zeros(50), np.full(50, 5.0)])
result = detect_mean(data)
print(result.cp_set)
```

The public API includes mean, variance, mean/variance, exponential, VAR,
linear, lasso, binomial, Poisson, quantile, GARCH, AR, ARMA, and ARIMA change
detection, plus rank and kernel transforms.

`detect_var(data, order=p)` accepts the raw multivariate time series and
constructs its lagged VAR design internally, matching the R interface. For an
already-constructed response/predictor matrix, use
`detect(data, family="mgaussian", p_response=q)` directly. The older
`var(data, order=p, p_response=q)` form for a pre-constructed matrix remains
accepted for compatibility.

`detect_arima(data, order=(p, d, q))` differences every candidate segment
independently in the shared native R/Python implementation. Returned change
points therefore use the original-series indices, and no cross-boundary
difference contaminates either adjacent segment. The likelihood is zero-mean
(`include_mean=False`), and `d=0` is identical to `detect_arma()`.

## Native build

Python packaging uses `scikit-build-core` and CMake. The source distribution
contains the shared C++ source and headers, so it builds without R or Bazel.
CMake fetches the pinned Armadillo headers and Abseil release. Linux builds
require BLAS/LAPACK development libraries (OpenBLAS is recommended), macOS
uses the system Accelerate framework, and Windows builds bundle the pinned
OpenBLAS DLL.

On Ubuntu/Debian, install the native prerequisites with:

```shell
sudo apt-get install g++ liblapack-dev libopenblas-dev
```

Then build and test an installed wheel:

```shell
python -m pip install -r requirements_lock.txt
python -m build --wheel
python -m pip install dist/fastcpd-*.whl
python -m pytest tests/test_fastcpd.py -m "not long"
```

Documentation: <https://x2r.io/fastcpd/python/>

Issues: <https://github.com/doccstat/fastcpd-py/issues>
