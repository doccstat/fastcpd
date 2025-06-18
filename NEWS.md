# fastcpd 0.18.0

*   Release wheel for Python on Linux and macOS.

# fastcpd 0.17.0

*   Introduce Python version.

# fastcpd 0.16.3

*   Move `fastcpd` to the top in the README.

# fastcpd 0.16.2

*   Remove natural language support.
*   Update mean variance cost function.
*   Use larger default `trim` value.

# fastcpd 0.16.1

*   Skip tests on CRAN.

# fastcpd 0.16.0

*   Performance improvement.
*   Less package dependencies.

# fastcpd 0.15.0

*   Remove constants.
*   Add parameter `variance_estimation` allowing users to specify the variance
    or covariance matrix if known.

# fastcpd 0.14.7

*   Simplify methods initialization.

# fastcpd 0.14.6

*   Update package comparison and remove packages breaking seed setting.

# fastcpd 0.14.5

*   Add pre-compiled vignettes.

# fastcpd 0.14.4

*   Simplify code.
*   Expose `fastcpd_impl` API for use in other packages.

# fastcpd 0.14.3

*   Skip certain tests on CRAN.

# fastcpd 0.14.2

*   Fix CRAN noSuggests error.

# fastcpd 0.14.1

*   Fix possible numeric issues on CRAN breaking updates of other packages.

# fastcpd 0.14.0

*   Correct mBIC.
*   Add citation.

# fastcpd 0.13.2

*   Remove package check in examples.
*   Use proper pruning coefficients for MBIC and MDL.
*   Make residuals matrices.
*   Use preprocess for faster mean change detection.
*   Update examples demonstrating faster mean change.
*   Replace comparison packages in vignettes with strange dependencies with
    results and `eval = FALSE`.

# fastcpd 0.13.1

*   Default no pruning for `lasso`.
*   Comment out `gfpop` due to https://github.com/doccstat/fastcpd/issues/10.

# fastcpd 0.13.0

*   Update all documentations.

# fastcpd 0.12.4

*   Customizable and pretty plots.

# fastcpd 0.12.3

*   Remove `pruning` parameter and replace with `convexity_coef = -Inf`.

# fastcpd 0.12.2

*   Update vignettes.

# fastcpd 0.12.1

*   Remove useless C++ codes.
*   Add more debug points in C++.
*   Add more examples for data `well_log`.
*   Add detection comparison for `well_log` data.
*   Add a variance estimator for median change.
*   Deprecate `winsorize_minval` and `winsorize_maxval`.

# fastcpd 0.12.0

*   Add Rice estimation for ARMA model variance estimation.
*   Add time comparison using Well-log data in vignettes.

# fastcpd 0.11.3

*   Add Rice estimator for mean change variance estimation.

# fastcpd 0.11.2

*   Export variance estimator function for linear models.

# fastcpd 0.11.1

*   Add package comparison with `CptNonPar`, `gfpop`, `InspectChangepoint`,
    `jointseg`, `Rbeast` and `VARDetect`.

# fastcpd 0.11.0

*   **Note**: From now on, MBIC is used as the default penalty selection for
    `beta` parameter.
*   Add penalty selection criteria using

    1. BIC: `(p + 1) * log(nrow(data)) / 2`
    1. Modified BIC: `(p + 2) * log(nrow(data)) / 2` with adjusted cost
       function.
    1. MDL: `(p + 2) * log(nrow(data)) / 2` with adjusted cost function.

    In the mean time, a numeric value can be passed to `beta` as well to
    explicitly specify the penalty for BIC.

# fastcpd 0.10.3

*   Add package check in examples and tests.

# fastcpd 0.10.2

*   Remove `bcp` according to

        Package ‘bcp’ was removed from the CRAN repository.

        Formerly available versions can be obtained from the archive.

        Archived on 2024-01-12 as email to the maintainer is undeliverable.

        A summary of the most recent check results can be obtained from the check results archive.

        Please use the canonical form https://CRAN.R-project.org/package=bcp to link to this page.

# fastcpd 0.10.1

*   Increase test coverage.
*   Use `interactive()` to check if the current R session is interactive.

# fastcpd 0.10.0

*   Add package comparison with other packages.

# fastcpd 0.9.9

*   Add small shiny app.

# fastcpd 0.9.8

*   Preliminary support for ARMA(p, q) model with parameter `order = c(p, q)`
    and family `"arma"`.
*   Add `fastcpd.arma` / `fastcpd_arma` for ARMA(p, q) model.
*   Add adaptive increasing `beta` values.

# fastcpd 0.9.7

*   Add `lower` and `upper` parameters to denote the lower and upper bounds of
    the parameters.
*   Add line search.
*   Add hardcoded ARMA(3, 2).

# fastcpd 0.9.6

*   Add `bitcoin` and `well_log` data.
*   Add residual calculation for mean family.
*   Add plots for bitcoin data.
*   Fix residual calculation for time series data when the seegments are
    too small.
*   Handle variance estimation errors.

# fastcpd 0.9.5

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

# fastcpd 0.9.4

*   Fix sanity check.
*   Add small AR(1) data in gallery.
*   Fix VAR(p) model bug.
*   Add VAR(2) example in Gallery.
*   Remove commented code.

# fastcpd 0.9.3

*   Deprecate "vanilla" family by `vanilla_percentage` parameter.
*   Add check utility functions.
*   Add MA(4) example.
*   Fix the bug when `beta` is updated but the old `beta` is still in use.
*   Remove tests estimating the variance in the "gaussian" family dynamically.
*   Merge `beta` updating into `get_segment_statistics`.

# fastcpd 0.9.2

*   Add gallery to vignettes.
*   Remove cheatsheets pdf from the package.
*   Use `forecast` package for ARIMA model.
*   Use `fGarch` package for GARCH model.

# fastcpd 0.9.1

*   Wrap `&&` around `||` by parentheses.
*   Add ma(4) example using custom cost function.
*   Add full support for AR(p), MA(q) and ARIMA(p, d, q) models.

# fastcpd 0.9.0

*   Submit for CRAN update.

# fastcpd 0.8.4

*   Add cheatsheets.
*   Refactor code and utilize the `cost_function_wrapper`.
*   Optimize warm start.

# fastcpd 0.8.3

*   Add `fastcpd.ts` / `fastcpd_ts` for time series data.
*   Fix pre segmengation bug for `lasso`.
*   Fix bug related to `vanilla_percentage` parameter for `lasso`.
*   Add tests with invalid family for `fastcpd.ts`.
*   Remove the `cp_only = TRUE` default when the family is "custom".
*   Improved plotting for "ar" and "var" families.
*   Add test coverage for `cp_only = TRUE` and `fastcpd_ts`.
*   Increase test coverage.
*   Provide user selection when `ggplot2` is not installed.

# fastcpd 0.8.2

*   Add cheatsheets WIP.
*   Add smaller examples test for penalized linear regression.

# fastcpd 0.8.1

*   Add new "ar" family for autoregressive models.
*   Add new "var" family for vector autoregressive models.

# fastcpd 0.8.0

*   Deal with the following:

        Due to the excessive calls to `glmnet` between R and C++,
        it is better to use the R implementation of `fastcpd` for lasso.

*   Separate the use of internal C++ cost functions and user-defined R cost
    functions.
*   Add Codecov Icicle plot in README.
*   Remove `cost_optim` and `cost_update` from `RcppExports.R`.
*   Estimate the variance in the "gaussian" family dynamically.

# fastcpd 0.7.6

*   Move default cost functions definition inside the `fastcpd` definition.
*   Define constant unordered set to store the family sets.
*   Avoid using `length(formals(cost))` to check the number of arguments of
    `cost` function.
*   Introduce an internal family "vanilla".

# fastcpd 0.7.5

*   Add variance estimation example in linear regression.
*   Update reference page.
*   Add validation for `family`.
*   Add user selection when `ggplot2` is not installed.
*   Add AR(1) using `forecast` example in the tests.

# fastcpd 0.7.4

*   Update website UI.
*   Update `fastcpd` documentation.

# fastcpd 0.7.3

*   Allow multiple response variables in the `formula`.
*   Add fastcpd logo to README.

# fastcpd 0.7.2

*   Add suggested package checking in tests.
*   Try to solve the *amazing* clang-ASAN error on CRAN:

        Error in dyn.load(file, DLLpath = DLLpath, ...) :
          unable to load shared object '/data/gannet/ripley/R/test-clang/mvtnorm/libs/mvtnorm.so':
          /data/gannet/ripley/R/test-clang/mvtnorm/libs/mvtnorm.so: undefined symbol: _ZNK7Fortran7runtime10Terminator5CrashEPKcz
        Calls: <Anonymous> ... asNamespace -> loadNamespace -> library.dynam -> dyn.load

# fastcpd 0.7.1

*   Add package citation.

# fastcpd 0.7.0

*   Remove C++ unit tests using catch and commented out the code since the new
    version of development version of Rcpp is not yet available on CRAN.
    Related pull request: https://github.com/RcppCore/Rcpp/pull/1274.
*   Add more documentation for `fastcpd` method.

# fastcpd 0.6.5

*   Add more experiments.

# fastcpd 0.6.4

*   Check warning messages in tests.
*   Further encapsulation of all FastcpdParameters members.

# fastcpd 0.6.3

*   Add CRAN release badge.

# fastcpd 0.6.2

*   Address CRAN comments.
*   Add more experiments.

# fastcpd 0.6.1

*   Address CRAN comments.

# fastcpd 0.6.0

*   Submit for CRAN release.

# fastcpd 0.5.7

*   Fix loss function for custom mean or variance change.

# fastcpd 0.5.6

*   Add stargazers in README.

# fastcpd 0.5.5

*   Add example and test for multivariate mean shift.
*   Add example and test for multivariate variance change.
*   Add example and test for multivariate mean and variance change.
*   Add test for linear regression with multi-dimensional responses.

# fastcpd 0.5.4

*   Fix bug when no change point is detected.

# fastcpd 0.5.3

*   Add more experiments but commented out for the sake of test time without
    affecting the test coverage.
*   Add examples in README.
*   Add CRAN manual using
    `R CMD Rd2pdf . --output=man/figures/manual.pdf --force --no-preview` from
    [stackoverflow](https://stackoverflow.com/a/50134068).
*   Add example for multiple epochs using custom cost functions.
*   Add table of contents to README.

# fastcpd 0.5.2

*   Add one-dimensional linear regression example with plot.

# fastcpd 0.5.1

*   Prepare for CRAN release.

# fastcpd 0.5.0

*   Rewrite the whole package in C++ except LASSO due to the excessive calls
    between R and C++ of `glmnet`.

# fastcpd 0.4.0

*   Add the transition from vanilla PELT to SeN by using `vanilla_percentage`
    parameter.

# fastcpd 0.3.3

*   Merge the implementation of vanilla PELT and SeN.
*   Encapsulate the implementation of binding new coefficients into the previous
    coefficients.
*   Rewrite `fastcpd` parameters updating in C++.

# fastcpd 0.3.2

*   Integrate initialization and update for `theta_hat`, `theta_sum` and
    `hessian`.
*   Combine theta estimation into a single function.
*   Add a parameter `vanilla_percentage` to denote the method switching between
    vanilla PETL and SeN.
*   Add more documentation to `cp_only` parameter.
*   Add preparation for merging vanilla PELT with SeN.

# fastcpd 0.3.1

*   Add examples as tests for `fastcpd`.
*   Rearrange C++ functions.
*   Add more precondition check.

# fastcpd 0.3.0

*   Bump test coverage for class methods of `fastcpd`.

# fastcpd 0.2.9

*   Fix Poisson regression bug related to `lfactorial`.

# fastcpd 0.2.8

*   Make penalized linear regression estimated coefficients output sparse.

# fastcpd 0.2.7

*   Fix mean change example bug.
*   Update documentation to redirect README to `pkgdown` generated webpage.
*   Add contact methods and ways to file a ticket.

# fastcpd 0.2.6

*   Add C++ sanity check for Logistic regression data, i.e. binomial family.
*   Add examples as tests for `fastcpd`.
*   Rename C++ source files to follow Unix convention.
*   Update documentation link in README.

# fastcpd 0.2.5

*   Hide internal functions from the documentation.
*   Export `fastcpd` class.

# fastcpd 0.2.4

*   Add column name for `thetas` slot in `fastcpd` class.
*   Fix plot where residuals and responses appear in the same plot.
*   Default `cp_only` to `FALSE`.
*   Remove residuals from `summary` method.

# fastcpd 0.2.3

*   Add missing examples for linear regression and LASSO.

# fastcpd 0.2.2

*   Add more examples to illustrate the use of the `fastcpd` function.
*   Indicating internal functions so that the users should not use them.

# fastcpd 0.2.1

*   Add more examples to the README.

# fastcpd 0.2.0

*   Added a `NEWS.md` file to track changes to the package.
