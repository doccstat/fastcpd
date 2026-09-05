
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Fast Change Point Detection <a href="https://fastcpd.xingchi.li"><img src="man/figures/logo.png" align="right" height="138" /></a>

[![Codecov test
coverage](https://codecov.io/gh/doccstat/fastcpd-r/branch/main/graph/badge.svg)](https://app.codecov.io/gh/doccstat/fastcpd-r?branch=main)
[![CodeFactor](https://www.codefactor.io/repository/github/doccstat/fastcpd-r/badge)](https://www.codefactor.io/repository/github/doccstat/fastcpd-r)
[![CRAN
status](https://www.r-pkg.org/badges/version-last-release/fastcpd)](https://cran.r-project.org/package=fastcpd)
[![doi](https://img.shields.io/badge/doi-10.48550/arXiv.2404.05933-green.svg)](https://doi.org/10.48550/arXiv.2404.05933)
[![R CMD check](https://github.com/doccstat/fastcpd-r/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/doccstat/fastcpd-r/actions/workflows/check-standard.yaml)
[![r-universe](https://doccstat.r-universe.dev/badges/fastcpd)](https://doccstat.r-universe.dev)
[![Python
version](https://img.shields.io/pypi/pyversions/fastcpd)](https://pypi.org/project/fastcpd/)
[![Python
package](https://img.shields.io/pypi/v/fastcpd)](https://pypi.org/project/fastcpd/)

## Documentation: [x2r.io](https://x2r.io/fastcpd/reference/)

Python and standalone C++ sources are published separately in
[`fastcpd-py`](https://github.com/doccstat/fastcpd-py) and
[`fastcpd-cpp`](https://github.com/doccstat/fastcpd-cpp).

Version 1.3.0 is the first coordinated R, Python, and standalone C++ source
interface. Portable built-in detectors share the same native algorithms,
defaults, seeded stochastic streams, change points, costs, parameters,
residual layout, and supported confidence diagnostics. Formula/data-frame
evaluation remains R-only. R closures/external pointers and C++
`std::function` callbacks are language-native extension points, while Python
rejects callbacks to preserve GIL-free detector execution. NumPy generator
streams and each language's result container also remain idiomatic extensions.
These boundaries keep the detector's compiled PELT and SEN paths free of
cross-language callback overhead.

## Three-language interface contract

| Capability | R | Python | Standalone C++ |
|---|---|---|
| Portable built-in detectors | `detect*()` on vectors, matrices, formulas, or data frames | `detect*()` on numeric array-like inputs | `detect*()` on Armadillo vectors or matrices |
| Deterministic numerical contract | Shared change points, costs, residuals, parameters, seeded KCP, and seeded bootstrap | Same | Same |
| Confidence intervals | `confint()` returns an R data frame | `confint()` returns dictionaries | `confint()` returns `ConfidenceInterval` values |
| Result ownership | S4 result retains its call and input | Frozen result retains read-only input and fit options | Lightweight result does not copy the input or fit options; pass them explicitly to `confint()` |
| Callback extensions | R closures or compiled external pointers; callable epoch schedules | Intentionally unavailable | Native `std::function` costs and epoch schedules |
| Compatibility aliases | R `fastcpd*()` spellings retained | Python short, `detect_*`, and `fastcpd_*` spellings retained | Canonical `detect_*` plus non-conflicting short names; `arma` and `garch` remain namespaces |

The standalone C++ source can build its dependencies for an in-tree build:

``` shell
git clone https://github.com/doccstat/fastcpd-cpp.git
cmake -S fastcpd-cpp -B fastcpd-cpp/build \
  -DFASTCPD_FETCH_DEPENDENCIES=ON \
  -DFASTCPD_INSTALL_CPP=OFF
cmake --build fastcpd-cpp/build --parallel
```

Installed CMake consumers use `find_package(fastcpd CONFIG REQUIRED)` and link
`fastcpd::fastcpd`. The installed package deliberately expects Armadillo (or
BLAS/LAPACK for a no-wrapper build) and Abseil 20260526 to be installed
separately; `FASTCPD_FETCH_DEPENDENCIES=ON` is a source-build convenience, not
a dependency-bundling install mode.

``` cpp
#include <fastcpd/fastcpd.h>

arma::mat data(100, 1, arma::fill::randn);
data.rows(50, 99) += 4.0;
fastcpd::Options fit_options;
fit_options.beta = 5.0;
fit_options.cost_adjustment = "BIC";
auto result = fastcpd::detect_mean(data, fit_options);

fastcpd::ConfidenceOptions interval_options;
interval_options.method = "profile";
interval_options.detector_options = fit_options;
auto intervals = fastcpd::confint(result, data, interval_options);
```

<details close>
<summary>
Installation
</summary>

``` r
# install.packages("pak")
pak::pak("doccstat/fastcpd-r")
# or install from CRAN
install.packages("fastcpd")
```

``` shell
# Install the Python package from PyPI
pip install fastcpd
```

</details>

## Comparison

``` r
set.seed(1)
n <- 10^8
mean_data <- c(rnorm(n / 2, 0, 1), rnorm(n / 2, 50, 1))
print(run_isolated(fastcpd::detect_mean(mean_data, cp_only = TRUE, variance_estimation = 1)))
#>    user  system elapsed 
#>   6.938   5.040  11.779
print(run_isolated(mosum::mosum(c(mean_data), G = 40)))
#>    user  system elapsed 
#>   9.034   6.854  15.964
print(run_isolated(changepoint::cpt.mean(mean_data, method = "PELT")))
#>    user  system elapsed 
#>  31.515   6.507  38.127
print(run_isolated(fpop::Fpop(mean_data, 2 * log(n))))
#>    user  system elapsed 
#>  44.527   2.978  47.632
```

![](man/figures/README-time-comparison-fastbench-1.png)<!-- -->

## References

- [fastcpd: Fast Change Point Detection in
  R](https://doi.org/10.48550/arXiv.2404.05933)
- [Sequential Gradient Descent and Quasi-Newton’s Method for
  Change-Point
  Analysis](https://proceedings.mlr.press/v206/zhang23b.html)

## FAQ

<details close>
<summary>
I countered problems related to gfortran on Mac OSX or Linux!
</summary>

The package should be able to install on Mac and any Linux distribution
without any problems if all the dependencies are installed. However, if
you encountered problems related to gfortran, it might be because
`RcppArmadillo` is not installed previously. Try [Mac OSX stackoverflow
solution](https://stackoverflow.com/a/72997915) or [Linux stackover
solution](https://stackoverflow.com/a/15540919) if you have trouble
installing `RcppArmadillo`.

</details>
<details close>
<summary>
We welcome contributions from everyone. Please follow the instructions
below to make contributions.
</summary>

1.  Fork the repo.

2.  Create a new branch from `main` branch.

3.  Make changes and commit them.

    1.  Please follow the [Google’s R style
        guide](https://google.github.io/styleguide/Rguide.html) for
        naming variables and functions.
    2.  If you are adding a new family of models with new cost functions
        with corresponding gradient and Hessian, please add them to
        `src/fastcpd_class_cost.cc` with proper example and tests in
        `vignettes/gallery.Rmd` and `tests/testthat/test-gallery.R`.
    3.  Add the family name to `src/fastcpd_constants.h`.
    4.  \[Recommended\] Add a new wrapper function in
        `R/fastcpd_wrappers.R` for the new family of models and move the
        examples to the new wrapper function as roxygen examples.
    5.  Add the new wrapper function to the corresponding section in
        `_pkgdown.yml`.

4.  Push the changes to your fork.

5.  Create a pull request.

6.  Make sure the pull request does not create new warnings or errors in
    `devtools::check()`.

</details>
<details close>
<summary>
Trouble installing Python package.
</summary>

Python headers are required to install the Python package. If you are
using Ubuntu, you can install the headers with:

``` shell
sudo apt install python3-dev
```

</details>
<details close>
<summary>
Encountered a bug or unintended behavior?
</summary>

1.  File a ticket at [GitHub
    Issues](https://github.com/doccstat/fastcpd-r/issues).
2.  Contact the authors specified in
    [DESCRIPTION](https://github.com/doccstat/fastcpd-r/blob/main/DESCRIPTION#L5-L10).

</details>
