# fastcpd 0.3.3

*   Merge the implementation of vanilla PELT and SeN.

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
