
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
#>   0.729   0.201   0.930
print(run_isolated(mosum::mosum(c(mean_data), G = 40)))
#>    user  system elapsed 
#>   1.219   0.667   1.893
print(run_isolated(changepoint::cpt.mean(mean_data, method = "PELT")))
#>    user  system elapsed 
#>   3.326   0.598   3.929
print(run_isolated(fpop::Fpop(mean_data, 2 * log(n))))
#>    user  system elapsed 
#>   3.953   0.262   4.218
```

![](https://raw.githubusercontent.com/doccstat/fastcpd-r/main/man/figures/README-time-comparison-fastbench-1.png)<!-- -->

### Python

``` python
import time

import fastcpd
import numpy as np
import ruptures
import sdt
import sdt.changepoint
import skchange
import skchange.detectors
import skchange.interval_scorers

rng = np.random.default_rng(1)
n = int(1e6)
x = np.r_[rng.normal(0, 1, n // 2), rng.normal(50, 1, n // 2)]
step = 10_000

start = time.perf_counter()
fastcpd.detect_mean(x, variance_estimation=1, cp_only=True)
print(f"fastcpd: {time.perf_counter() - start:.3f} s")

start = time.perf_counter()
skchange.detectors.PELT(
    cost=skchange.interval_scorers.L2Cost(), penalty=2 * np.log(n), step_size=step
).fit_predict(x.reshape(-1, 1))
print(f"skchange: {time.perf_counter() - start:.3f} s")

start = time.perf_counter()
sdt.changepoint.Pelt(cost="l2", min_size=step, jump=step).find_changepoints(
    x, penalty=2 * np.log(n)
)
print(f"sdt-python: {time.perf_counter() - start:.3f} s")

start = time.perf_counter()
ruptures.Pelt(model="l2", min_size=step, jump=step).fit(x).predict(
    pen=2 * np.log(n)
)
print(f"ruptures: {time.perf_counter() - start:.3f} s")
```

    #> fastcpd: 0.061 s
    #> skchange: 0.316 s
    #> sdt-python: 5.989 s
    #> ruptures: 0.852 s

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
