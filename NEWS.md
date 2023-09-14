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
