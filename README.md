
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Fast Change Point Detection <a href="https://fastcpd.xingchi.li"><img src="man/figures/logo.png" align="right" height="138" /></a>

[![Codecov test
coverage](https://codecov.io/gh/doccstat/fastcpd/branch/main/graph/badge.svg)](https://app.codecov.io/gh/doccstat/fastcpd?branch=main)
[![CodeFactor](https://www.codefactor.io/repository/github/doccstat/fastcpd/badge)](https://www.codefactor.io/repository/github/doccstat/fastcpd)
[![CRAN
status](https://www.r-pkg.org/badges/version-last-release/fastcpd)](https://cran.r-project.org/package=fastcpd)
[![doi](https://img.shields.io/badge/doi-10.48550/arXiv.2404.05933-green.svg)](https://doi.org/10.48550/arXiv.2404.05933)
[![R-CMD-check.yaml](https://github.com/doccstat/fastcpd/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/doccstat/fastcpd/actions)
[![r-universe](https://doccstat.r-universe.dev/badges/fastcpd)](https://doccstat.r-universe.dev)
[![Python
version](https://img.shields.io/pypi/pyversions/fastcpd)](https://pypi.org/project/fastcpd/)
[![Python
package](https://img.shields.io/pypi/v/fastcpd)](https://pypi.org/project/fastcpd/)

## Documentation

- R documentation: [fastcpd.xingchi.li](https://fastcpd.xingchi.li)
- Python documentation:
  [fastcpd.xingchi.li/python](https://fastcpd.xingchi.li/python/)

## Installation

### R

``` r
# install.packages("pak")
pak::pak("doccstat/fastcpd")
# or install from CRAN
install.packages("fastcpd")
```

### Python WIP

``` shell
# python -m ensurepip --upgrade
pip install .
# or install from PyPI
pip install fastcpd
```

## Usage

### R

``` r
set.seed(1)
n <- 1000
x <- rep(0, n + 3)
for (i in 1:600) {
  x[i + 3] <- 0.6 * x[i + 2] - 0.2 * x[i + 1] + 0.1 * x[i] + rnorm(1, 0, 3)
}
for (i in 601:1000) {
  x[i + 3] <- 0.3 * x[i + 2] + 0.4 * x[i + 1] + 0.2 * x[i] + rnorm(1, 0, 3)
}
result <- fastcpd::fastcpd.ar(x[3 + seq_len(n)], 3, r.progress = FALSE)
summary(result)
#> 
#> Call:
#> fastcpd::fastcpd.ar(data = x[3 + seq_len(n)], order = 3, r.progress = FALSE)
#> 
#> Change points:
#> 614 
#> 
#> Cost values:
#> 2754.116 2038.945 
#> 
#> Parameters:
#>     segment 1 segment 2
#> 1  0.57120256 0.2371809
#> 2 -0.20985108 0.4031244
#> 3  0.08221978 0.2290323
plot(result)
```

![](man/figures/README-ar3-1.png)<!-- -->

### Python WIP

``` python
import fastcpd.segmentation
from numpy import concatenate
from numpy.random import normal, multivariate_normal
covariance_mat = [[100, 0, 0], [0, 100, 0], [0, 0, 100]]
data = concatenate((multivariate_normal([0, 0, 0], covariance_mat, 300),
                    multivariate_normal([50, 50, 50], covariance_mat, 400),
                    multivariate_normal([2, 2, 2], covariance_mat, 300)))
fastcpd.segmentation.mean(data)

import fastcpd.variance_estimation
fastcpd.variance_estimation.mean(data)
```

### Comparison

``` r
set.seed(1)
n <- 10^8
mean_data <- c(rnorm(n / 2, 0, 1), rnorm(n / 2, 50, 1))
run_isolated <- function(expr) {
  callr::r(function(e, n) {
    set.seed(1)
    mean_data <- c(rnorm(n / 2, 0, 1), rnorm(n / 2, 50, 1))
    system.time(eval(e))
  }, args = list(e = substitute(expr), n = n))
}
print(run_isolated(fastcpd::fastcpd.mean(mean_data, r.progress = FALSE, cp_only = TRUE, variance_estimation = 1)))
#>    user  system elapsed 
#>   8.593   6.494  14.851
print(run_isolated(mosum::mosum(c(mean_data), G = 40)))
#>    user  system elapsed 
#>   9.077   6.882  16.029
print(run_isolated(fpop::Fpop(mean_data, 2 * log(n))))
#>    user  system elapsed 
#>  44.990   3.088  48.226
print(run_isolated(changepoint::cpt.mean(mean_data, method = "PELT")))
#>    user  system elapsed 
#>  31.516   6.531  38.159
```

![](man/figures/README-time-comparison-fastbench-1.png)<!-- -->

## Cheatsheet

[![fastcpd
cheatsheet](man/figures/cheatsheets.png)](https://github.com/doccstat/fastcpd/blob/main/man/figures/cheatsheets.pdf)

### References

- [fastcpd: Fast Change Point Detection in
  R](https://doi.org/10.48550/arXiv.2404.05933)
- [Sequential Gradient Descent and Quasi-Newton’s Method for
  Change-Point
  Analysis](https://proceedings.mlr.press/v206/zhang23b.html)

## FAQ

<details close>
<summary>
Should I install suggested packages?
</summary>

The suggested packages are not required for the main functionality of
the package. They are only required for the vignettes. If you want to
learn more about the package comparison and other vignettes, you could
either check out vignettes on
[CRAN](https://CRAN.R-project.org/package=fastcpd) or [pkgdown generated
documentation](https://fastcpd.xingchi.li/articles/).

</details>
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
    Issues](https://github.com/doccstat/fastcpd/issues).
2.  Contact the authors specified in
    [DESCRIPTION](https://github.com/doccstat/fastcpd/blob/main/DESCRIPTION#L5-L10).

</details>

## Stargazers over time

[![Stargazers over
time](https://starchart.cc/doccstat/fastcpd.svg)](https://starchart.cc/doccstat/fastcpd)
