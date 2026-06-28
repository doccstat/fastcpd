# Custom logistic regression model

We try to reproduce the logistic regression models with custom cost
functions and show that the results are similar to the built-in logistic
regression models.

The built-in logistic regression model in
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) is
implemented with the help of the **fastglm** package. The **fastglm**
package utilizes the iteratively reweighted least squares with the
step-halving approach to help safeguard against convergence issues. If a
custom cost function is used with gradient descent, we should expect the
results will be similar to the built-in logistic regression model.

Specifying the `cost`, `cost_gradient` and `cost_hessian` parameters
below, we can obtain similar results as the built-in logistic regression
model.

``` r

set.seed(1)
x <- matrix(rnorm(2000, 0, 1), ncol = 5)
theta <- rbind(rnorm(5, 0, 1), rnorm(5, 2, 1))
y <- c(
  rbinom(250, 1, 1 / (1 + exp(-x[1:250, ] %*% theta[1, ]))),
  rbinom(150, 1, 1 / (1 + exp(-x[251:400, ] %*% theta[2, ])))
)
binomial_data <- data.frame(y = y, x = x)

result <- fastcpd.binomial(cbind(y, x), cost_adjustment = NULL)
summary(result)
#> 
#> Call:
#> fastcpd.binomial(data = cbind(y, x), cost_adjustment = NULL)
#> 
#> Change points:
#> 250 
#> 
#> Cost values:
#> 99.40994 36.89336 
#> 
#> Parameters:
#>    segment 1  segment 2
#> 1 -1.0096300 3.46131278
#> 2 -1.7740956 2.10636392
#> 3  1.2736663 0.74124627
#> 4  0.3955351 0.07297997
#> 5 -0.1047536 2.00553153
```

``` r

logistic_loss <- function(data, theta) {
  x <- data[, -1]
  y <- data[, 1]
  u <- x %*% theta
  nll <- -y * u + log(1 + exp(u))
  nll[u > 10] <- -y[u > 10] * u[u > 10] + u[u > 10]
  sum(nll)
}
logistic_gradient <- function(data, theta) {
  x <- data[nrow(data), -1]
  y <- data[nrow(data), 1]
  c(-(y - 1 / (1 + exp(-x %*% theta)))) * x
}
logistic_hessian <- function(data, theta) {
  x <- data[nrow(data), -1]
  prob <- 1 / (1 + exp(-x %*% theta))
  (x %o% x) * c((1 - prob) * prob)
}
result <- fastcpd(
  y ~ . - 1, binomial_data, epsilon = 1e-5, cost = logistic_loss,
  cost_gradient = logistic_gradient, cost_hessian = logistic_hessian
)
summary(result)
#> 
#> Call:
#> fastcpd(formula = y ~ . - 1, data = binomial_data, cost = logistic_loss, 
#>     cost_gradient = logistic_gradient, cost_hessian = logistic_hessian, 
#>     epsilon = 1e-05)
#> 
#> Change points:
#> 252 
#> 
#> Parameters:
#>   segment 1 segment 2
#> 1         0         0
#> 2         0         0
#> 3         0         0
#> 4         0         0
#> 5         0         0
```

Note that the result obtained through custom cost functions is inferior
compared to the one obtained through built-in models. We remark that the
results can be improved with extra parameters already provided in the
package. The detailed discussion of several advanced usages of the
package can be found in [Advanced
examples](https://fastcpd.xingchi.li/articles/examples-advanced.md).

## Compiled (C++) custom cost functions

The custom cost functions above are R closures: `cost`, `cost_gradient`
and `cost_hessian` are called from C++ once per candidate segment via
`Rcpp::Function`, which incurs R call overhead. For performance-critical
use cases,
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) also
accepts **pre-compiled C++ cost functions**, passed through the very
same `cost` / `cost_gradient` / `cost_hessian` arguments as an
`externalptr` built with `Rcpp::XPtr`. The compiled and R-closure paths
compute the same mathematical cost and a `family = "custom"` run accepts
either interchangeably –
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) detects
which kind of SEXP it received (`TYPEOF(x) == EXTPTRSXP`) and dispatches
accordingly, with no other change to the call.

To build a compiled cost, write a C++ function matching one of the two
expected signatures – a PELT-style `cost(data)`:

``` cpp
double my_cost_pelt(arma::mat const& data);
typedef double (*CostPeltFnPtr)(arma::mat const&);
```

or a SeGD-style `cost(data, theta)` (optionally paired with a gradient
and a Hessian):

``` cpp
double my_cost_sen(arma::mat const& data, arma::colvec const& theta);
typedef double (*CostSenFnPtr)(arma::mat const&, arma::colvec const&);

arma::colvec my_cost_gradient(arma::mat const& data, arma::colvec const& theta);
typedef arma::colvec (*CostGradientFnPtr)(arma::mat const&, arma::colvec const&);

arma::mat my_cost_hessian(arma::mat const& data, arma::colvec const& theta);
typedef arma::mat (*CostHessianFnPtr)(arma::mat const&, arma::colvec const&);
```

Wrap the function pointer in an `Rcpp::XPtr`, tagging it with the
literal string `fastcpd` expects for that role (`"fastcpd_cost_pelt"`,
`"fastcpd_cost_sen"`, `"fastcpd_cost_gradient"` or
`"fastcpd_cost_hessian"` – this tag is how `fastcpd` validates that the
pointer was built for the right purpose), and export a function that
returns it:

``` cpp
// [[Rcpp::export]]
SEXP make_my_cost_pelt_xptr() {
  Rcpp::XPtr<CostPeltFnPtr> xptr(
      new CostPeltFnPtr(&my_cost_pelt), true,
      Rcpp::wrap("fastcpd_cost_pelt"));
  return xptr;
}
```

Because external pointers carry no
[`formals()`](https://rdrr.io/r/base/formals.html),
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) cannot
infer whether a compiled `cost` is PELT-style (one argument) or
SeGD-style (two arguments) the way it does for R closures. Instead, set
an integer `fastcpd_cost_arity` attribute on the returned pointer – `1L`
for `cost(data)`, `2L` for `cost(data, theta)`:

``` r

my_cost_pelt_xptr <- make_my_cost_pelt_xptr()
attr(my_cost_pelt_xptr, "fastcpd_cost_arity") <- 1L

result <- fastcpd(
  ~ . - 1, data.frame(x = data), family = "custom", cost = my_cost_pelt_xptr
)
summary(result)
```

A compiled `cost` can be used on its own (covering both the PELT-style
and SeGD-style cases above), or paired with compiled `cost_gradient` and
`cost_hessian` `externalptr`s the same way – including the
gradient/Hessian- driven PELT warm start (`fastcpd_cost_arity = 1L`),
which is now solved by a from-scratch Brent / bound-constrained BFGS
optimiser working directly through the unified
`cost`/`cost_gradient`/`cost_hessian` wrappers, with no R-level callback
and therefore no restriction on whether `cost` is an R closure or a
compiled `externalptr`.

## Notes

This document is generated by the following code:

``` shell
R -e 'knitr::knit("vignettes/examples-custom-model.Rmd.original", output = "vignettes/examples-custom-model.Rmd")'
```

## Appendix: all code snippets

``` r

knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>", eval = TRUE, warning = FALSE
)
library(fastcpd)
set.seed(1)
x <- matrix(rnorm(2000, 0, 1), ncol = 5)
theta <- rbind(rnorm(5, 0, 1), rnorm(5, 2, 1))
y <- c(
  rbinom(250, 1, 1 / (1 + exp(-x[1:250, ] %*% theta[1, ]))),
  rbinom(150, 1, 1 / (1 + exp(-x[251:400, ] %*% theta[2, ])))
)
binomial_data <- data.frame(y = y, x = x)

result <- fastcpd.binomial(cbind(y, x), cost_adjustment = NULL)
summary(result)
logistic_loss <- function(data, theta) {
  x <- data[, -1]
  y <- data[, 1]
  u <- x %*% theta
  nll <- -y * u + log(1 + exp(u))
  nll[u > 10] <- -y[u > 10] * u[u > 10] + u[u > 10]
  sum(nll)
}
logistic_gradient <- function(data, theta) {
  x <- data[nrow(data), -1]
  y <- data[nrow(data), 1]
  c(-(y - 1 / (1 + exp(-x %*% theta)))) * x
}
logistic_hessian <- function(data, theta) {
  x <- data[nrow(data), -1]
  prob <- 1 / (1 + exp(-x %*% theta))
  (x %o% x) * c((1 - prob) * prob)
}
result <- fastcpd(
  y ~ . - 1, binomial_data, epsilon = 1e-5, cost = logistic_loss,
  cost_gradient = logistic_gradient, cost_hessian = logistic_hessian
)
summary(result)
my_cost_pelt_xptr <- make_my_cost_pelt_xptr()
attr(my_cost_pelt_xptr, "fastcpd_cost_arity") <- 1L

result <- fastcpd(
  ~ . - 1, data.frame(x = data), family = "custom", cost = my_cost_pelt_xptr
)
summary(result)
```
