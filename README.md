
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Fast Change Point Detection <a href="https://fastcpd.xingchi.li"><img src="https://raw.githubusercontent.com/doccstat/fastcpd-r/main/man/figures/logo.png" align="right" height="138" /></a>

[![Codecov test
coverage](https://codecov.io/gh/doccstat/fastcpd-r/branch/main/graph/badge.svg)](https://app.codecov.io/gh/doccstat/fastcpd-r?branch=main)
[![CodeFactor](https://www.codefactor.io/repository/github/doccstat/fastcpd-r/badge)](https://www.codefactor.io/repository/github/doccstat/fastcpd-r)
[![CRAN
status](https://www.r-pkg.org/badges/version-last-release/fastcpd)](https://cran.r-project.org/package=fastcpd)
[![doi](https://img.shields.io/badge/doi-10.48550/arXiv.2404.05933-green.svg)](https://doi.org/10.48550/arXiv.2404.05933)
[![R CMD
check](https://github.com/doccstat/fastcpd-r/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/doccstat/fastcpd-r/actions/workflows/check-standard.yaml)
[![r-universe](https://doccstat.r-universe.dev/badges/fastcpd)](https://doccstat.r-universe.dev)
[![Python
version](https://img.shields.io/pypi/pyversions/fastcpd)](https://pypi.org/project/fastcpd/)
[![Python
package](https://img.shields.io/pypi/v/fastcpd)](https://pypi.org/project/fastcpd/)

## Documentation: [x2r.io](https://x2r.io/fastcpd/reference/)

Python and standalone C++ sources are published separately in
[`fastcpd-py`](https://github.com/doccstat/fastcpd-py) and
[`fastcpd-cpp`](https://github.com/doccstat/fastcpd-cpp).

<details>
<summary>
Installation: R, Python, and C++
</summary>

R package:

``` r
install.packages("fastcpd")
```

Python package:

``` shell
python -m pip install fastcpd
```

C++ library (requires Armadillo and Abseil 20260526 or newer):

``` shell
git clone https://github.com/doccstat/fastcpd-cpp.git
cmake -S fastcpd-cpp -B fastcpd-cpp/build -DFASTCPD_BUILD_EXAMPLES=OFF
cmake --build fastcpd-cpp/build --parallel
cmake --install fastcpd-cpp/build --prefix fastcpd-install
```

</details>

## Comparison

### R

``` r
set.seed(1)
n <- 10^7
mean_data <- c(rnorm(n / 2, 0, 1), rnorm(n / 2, 50, 1))
print(run_isolated(fastcpd::detect_mean(mean_data, cp_only = TRUE, variance_estimation = 1)))
#>    user  system elapsed 
#>   0.745   0.194   0.981
print(run_isolated(mosum::mosum(c(mean_data), G = 40)))
#>    user  system elapsed 
#>   1.282   0.696   2.236
print(run_isolated(changepoint::cpt.mean(mean_data, method = "PELT")))
#>    user  system elapsed 
#>   3.381   0.665   4.418
print(run_isolated(fpop::Fpop(mean_data, 2 * log(n))))
#>    user  system elapsed 
#>   4.027   0.299   4.483
```

![](https://raw.githubusercontent.com/doccstat/fastcpd-r/main/man/figures/README-time-comparison-fastbench-1.png)<!-- -->

### Python

``` python
import numpy as np
import ruptures as rpt
from fastcpd import detect_mean
from sdt.changepoint import Pelt as SdtPelt
from skchange.detectors import PELT
from skchange.interval_scorers import L2Cost

rng = np.random.default_rng(1)
n = 1_000
x = np.r_[rng.normal(0, 1, n // 2), rng.normal(50, 1, n // 2)]


benchmark("fastcpd", lambda: detect_mean(x, variance_estimation=1, cp_only=True))
benchmark("skchange", lambda: PELT(cost=L2Cost(), penalty=2 * np.log(n), min_segment_length=1).fit_predict(x.reshape(-1, 1)))
benchmark("sdt-python", lambda: SdtPelt(cost="l2", min_size=1, jump=1).find_changepoints(x, penalty=2 * np.log(n)))
benchmark("ruptures", lambda: rpt.Pelt(model="l2", min_size=1, jump=1).fit(x).predict(pen=2 * np.log(n)))
```

![](https://raw.githubusercontent.com/doccstat/fastcpd-r/main/man/figures/README-time-comparison-python-plot-1.png)<!-- -->

### C++

Native fastcpd and fpop on Linux ARM64, with 1,000,000 observations.
[Source and build
instructions](https://github.com/doccstat/fastcpd-r/tree/main/tools/readme-comparison).

![](https://raw.githubusercontent.com/doccstat/fastcpd-r/main/man/figures/README-time-comparison-cpp-plot-1.png)<!-- -->

## References

- [fastcpd: Fast Change Point Detection in
  R](https://doi.org/10.48550/arXiv.2404.05933)
- [Sequential Gradient Descent and Quasi-Newton’s Method for
  Change-Point
  Analysis](https://proceedings.mlr.press/v206/zhang23b.html)
