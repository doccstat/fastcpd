# Find change points efficiently

`fastcpd()` takes in formulas, data, families and extra parameters and
returns a
[fastcpd](https://fastcpd.xingchi.li/reference/fastcpd-class.md) object.

## Usage

``` r
fastcpd(
  formula = y ~ . - 1,
  data,
  beta = "MBIC",
  cost_adjustment = "MBIC",
  family = NULL,
  cost = NULL,
  cost_gradient = NULL,
  cost_hessian = NULL,
  line_search = c(1),
  lower = rep(-Inf, p),
  upper = rep(Inf, p),
  pruning_coef = 0,
  segment_count = 10,
  trim = 0.05,
  momentum_coef = 0,
  multiple_epochs = function(x) 0,
  epsilon = 1e-10,
  order = c(0, 0, 0),
  p = ncol(data) - 1,
  variance_estimation = NULL,
  cp_only = FALSE,
  vanilla_percentage = 0,
  warm_start = FALSE,
  ...
)
```

## Arguments

- formula:

  A formula object specifying the model to be fitted. The (optional)
  response variable should be on the LHS of the formula, while the
  covariates should be on the RHS. The naming of variables used in the
  formula should be consistent with the column names in the data frame
  provided in `data`. The intercept term should be removed from the
  formula. The response variable is not needed for mean/variance change
  models and time series models. By default, an intercept column will be
  added to the data, similar to the
  [`lm()`](https://rdrr.io/r/stats/lm.html) function. Thus, it is
  suggested that users should remove the intercept term by appending
  `- 1` to the formula. Note that the
  [fastcpd.family](https://fastcpd.xingchi.li/reference/fastcpd_family.md)
  functions do not require a formula input.

- data:

  A data frame of dimension \\T \times d\\ containing the data to be
  segmented (where each row denotes a data point \\z_t \in
  \mathbb{R}^d\\ for \\t = 1, \ldots, T\\) is required in the main
  function, while a matrix or a vector input is also accepted in the
  [fastcpd.family](https://fastcpd.xingchi.li/reference/fastcpd_family.md)
  functions.

- beta:

  Penalty criterion for the number of change points. This parameter
  takes a string value of `"BIC"`, `"MBIC"`, `"MDL"` or a numeric value.
  If a numeric value is provided, the value will be used as the penalty.
  By default, the mBIC criterion is used, where \\\beta = (p + 2)
  \log(T) / 2\\. This parameter usage should be paired with
  `cost_adjustment` described below. Discussions about the penalty
  criterion can be found in the references.

- cost_adjustment:

  Cost adjustment criterion. It can be `"BIC"`, `"MBIC"`, `"MDL"` or
  `NULL`. By default, the cost adjustment criterion is set to be
  `"MBIC"`. The `"MBIC"` and `"MDL"` criteria modify the cost function
  by adding a negative adjustment term to the cost function. `"BIC"` or
  `NULL` does not modify the cost function. Details can in found in the
  references.

- family:

  Family class of the change point model. It can be `"mean"` for mean
  change, `"variance"` for variance change, `"meanvariance"` for mean
  and/or variance change, `"lm"` for linear regression, `"binomial"` for
  logistic regression, `"poisson"` for Poisson regression, `"lasso"` for
  penalized linear regression, `"ar"` for AR(\\p\\) models, `"arma"` for
  ARMA(\\p\\, \\q\\) models, `"arima"` for ARIMA(\\p\\, \\d\\, \\q\\)
  models, `"garch"` for GARCH(\\p\\, \\q\\) models, `"var"` for
  VAR(\\p\\) models and `"custom"` for user-specified custom models.
  Omitting this parameter is the same as specifying the parameter to be
  `"custom"` or `NULL`, in which case, users must specify the custom
  cost function.

- cost:

  Cost function to be used. `cost`, `cost_gradient`, and `cost_hessian`
  should not be specified at the same time with `family` as built-in
  families have cost functions implemented in C++ to provide better
  performance. If not specified, the default is the negative
  log-likelihood for the corresponding family. Custom cost functions can
  be provided in the following two formats:

  - `cost = function(data) {...}`

  - `cost = function(data, theta) {...}`

  Users can specify a loss function using the second format that will be
  used to calculate the cost value. In both formats, the input data is a
  subset of the original data frame in the form of a matrix (a matrix
  with a single column in the case of a univariate data set). In the
  first format, the specified cost function directly calculates the cost
  value. `fastcpd()` performs the vanilla PELT algorithm, and
  `cost_gradient` and `cost_hessian` should not be provided since no
  parameter updating is necessary for vanilla PELT. In the second
  format, the loss function \\\sum\_{i = s}^t l(z_i, \theta)\\ is
  provided, which has to be optimized over the parameter \\\theta\\ to
  obtain the cost value. A detailed discussion about the custom cost
  function usage can be found in the references.

- cost_gradient:

  Gradient of the custom cost function. Example usage:

      cost_gradient = function(data, theta) {
        ...
        return(gradient)
      }

  The gradient function takes two inputs, the first being a matrix
  representing a segment of the data, similar to the format used in the
  `cost` function, and the second being the parameter that needs to be
  optimized. The gradient function returns the value of the gradient of
  the loss function, i.e., \\\sum\_{i = s}^t \nabla l(z_i, \theta)\\.

- cost_hessian:

  Hessian of the custom loss function. The Hessian function takes two
  inputs, the first being a matrix representing a segment of the data,
  similar to the format used in the `cost` function, and the second
  being the parameter that needs to be optimized. The gradient function
  returns the Hessian of the loss function, i.e., \\\sum\_{i = s}^t
  \nabla^2 l(z_i, \theta)\\.

- line_search:

  If a vector of numeric values is provided, a line search will be
  performed to find the optimal step size for each update. Detailed
  usage of `line_search` can be found in the references.

- lower:

  Lower bound for the parameters. Used to specify the domain of the
  parameters after each gradient descent step. If not specified, the
  lower bound is set to be `-Inf` for all parameters. `lower` is
  especially useful when the estimated parameters take only positive
  values, such as the noise variance.

- upper:

  Upper bound for the parameters. Used to specify the domain of the
  parameters after each gradient descent step. If not specified, the
  upper bound is set to be `Inf` for all parameters.

- pruning_coef:

  Pruning coefficient \$c_0\$ used in the pruning step of the PELT
  algorithm with the default value 0. If `cost_adjustment` is specified
  as `"MBIC"`, an adjustment term \\p\log(2)\\ will be added to the
  pruning coefficient. If `cost_adjustment` is specified as `"MDL"`, an
  adjustment term \\p\log_2(2)\\ will be added to the pruning
  coefficient. Detailed discussion about the pruning coefficient can be
  found in the references.

- segment_count:

  An initial guess of the number of segments. If not specified, the
  initial guess of the number of segments is 10. The initial guess
  affects the initial estimates of the parameters in SeGD.

- trim:

  Trimming for the boundary change points so that a change point close
  to the boundary will not be counted as a change point. This parameter
  also specifies the minimum distance between two change points. If
  several change points have mutual distances smaller than
  `trim * nrow(data)`, those change points will be merged into one
  single change point. The value of this parameter should be between 0
  and 1.

- momentum_coef:

  Momentum coefficient to be applied to each update. This parameter is
  used when the loss function is bad-shaped so that maintaining a
  momentum from previous update is desired. Default value is 0, meaning
  the algorithm doesn't maintain a momentum by default.

- multiple_epochs:

  A function can be specified such that an adaptive number of multiple
  epochs can be utilized to improve the algorithm's performance.
  `multiple_epochs` is a function of the length of the data segment. The
  function returns an integer indicating how many epochs should be
  performed apart from the default update. By default, the function
  returns zero, meaning no multiple epochs will be used to update the
  parameters. Example usage:

      multiple_epochs = function(segment_length) {
        if (segment_length < 100) 1
        else 0
      }

  This function will let SeGD perform parameter updates with an
  additional epoch for each segment with a length less than 100 and no
  additional epoch for segments with lengths greater or equal to 100.

- epsilon:

  Epsilon to avoid numerical issues. Only used for the Hessian
  computation in Logistic Regression and Poisson Regression.

- order:

  Order of the AR(\\p\\), VAR(\\p\\) or ARIMA(\\p\\, \\d\\, \\q\\)
  model.

- p:

  Number of covariates in the model. If not specified, the number of
  covariates will be inferred from the data, i.e., `p = ncol(data) - 1`.
  This parameter is superseded by `order` in the case of time series
  models: "ar", "var", "arima".

- variance_estimation:

  An estimate of the variance / covariance matrix for the data. If not
  specified, the variance / covariance matrix will be estimated using
  the data.

- cp_only:

  If `TRUE`, only the change points are returned. Otherwise, the cost
  function values together with the estimated parameters for each
  segment are also returned. By default the value is set to be `FALSE`
  so that `plot` can be used to visualize the results for a built-in
  model. `cp_only` has some performance impact on the algorithm, since
  the cost values and estimated parameters for each segment need to be
  calculated and stored. If the users are only interested in the change
  points, setting `cp_only` to be `TRUE` will help with the
  computational cost.

- vanilla_percentage:

  The parameter \\v\\ is between zero and one. For each segment, when
  its length is no more than \\vT\\, the cost value will be computed by
  performing an exact minimization of the loss function over the
  parameter. When its length is greater than \\vT\\, the cost value is
  approximated through SeGD. Therefore, this parameter induces an
  algorithm that can be interpreted as an interpolation between dynamic
  programming with SeGD (\\v = 0\\) and the vanilla PELT (\\v = 1\\).
  The readers are referred to the references for more details.

- warm_start:

  If `TRUE`, the algorithm will use the estimated parameters from the
  previous segment as the initial value for the current segment. This
  parameter is only used for the `"glm"` families.

- ...:

  Other parameters for specific models.

  - `include.mean` is used to determine if a mean/intercept term should
    be included in the ARIMA(\\p\\, \\d\\, \\q\\) or GARCH(\\p\\, \\q\\)
    models.

  - `r.progress` is used to control the progress bar. By default the
    progress bar will be shown. To disable it, set `r.progress = FALSE`.

  - `p.response` is used to specify the number of response variables.
    This parameter is especially useful for linear models with
    multivariate responses.

## Value

A [fastcpd](https://fastcpd.xingchi.li/reference/fastcpd-class.md)
object.

## Gallery

<https://github.com/doccstat/fastcpd/tree/main/tests/testthat/examples>

## References

Xingchi Li, Xianyang Zhang (2024). “fastcpd: Fast Change Point Detection
in R.” *arXiv:2404.05933*, <https://arxiv.org/abs/2404.05933>.

Xianyang Zhang, Trisha Dawn (2023). “Sequential Gradient Descent and
Quasi-Newton's Method for Change-Point Analysis.” In Ruiz, Francisco,
Dy, Jennifer, van de Meent, Jan-Willem (eds.), *Proceedings of The 26th
International Conference on Artificial Intelligence and Statistics*,
volume 206 series Proceedings of Machine Learning Research, 1129-1143.

## See also

[fastcpd.family](https://fastcpd.xingchi.li/reference/fastcpd_family.md)
for the family-specific function;
[`plot.fastcpd()`](https://fastcpd.xingchi.li/reference/plot.md) for
plotting the results,
[`summary.fastcpd()`](https://fastcpd.xingchi.li/reference/summary.md)
for summarizing the results.

## Examples

``` r
if (requireNamespace("mvtnorm", quietly = TRUE)) {
  set.seed(1)
  n <- 200
  p <- 4
  d <- 2
  x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
  theta_1 <- matrix(runif(8, -3, -1), nrow = p)
  theta_2 <- matrix(runif(8, -1, 3), nrow = p)
  y <- rbind(
    x[1:125, ] %*% theta_1 + mvtnorm::rmvnorm(125, rep(0, d), 3 * diag(d)),
    x[126:n, ] %*% theta_2 + mvtnorm::rmvnorm(75, rep(0, d), 3 * diag(d))
  )
  result_mlm <- fastcpd(
    cbind(y.1, y.2) ~ . - 1, cbind.data.frame(y = y, x = x), family = "lm"
  )
  summary(result_mlm)
}
#> 
#> Call:
#> fastcpd(formula = cbind(y.1, y.2) ~ . - 1, data = cbind.data.frame(y = y, 
#>     x = x), family = "lm")
#> 
#> Change points:
#> 125 
#> 
#> Cost values:
#> 522.26 317.5959 
#> 
#> Parameters:
#>   segment 1  segment 2
#> 1 -2.870872  0.1521972
#> 2 -2.946976 -0.4559850
#> 3 -2.971442  0.6369559
#> 4 -1.257543  1.5949148
#> 5 -1.501633  0.3805461
#> 6 -2.737303  1.9325768
#> 7 -2.034249  2.2404699
#> 8 -2.438544  2.2219857
if (
  requireNamespace("mvtnorm", quietly = TRUE) &&
    requireNamespace("stats", quietly = TRUE)
) {
  set.seed(1)
  n <- 400 + 300 + 500
  p <- 5
  x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
  theta <- rbind(
    mvtnorm::rmvnorm(1, mean = rep(0, p - 3), sigma = diag(p - 3)),
    mvtnorm::rmvnorm(1, mean = rep(5, p - 3), sigma = diag(p - 3)),
    mvtnorm::rmvnorm(1, mean = rep(9, p - 3), sigma = diag(p - 3))
  )
  theta <- cbind(theta, matrix(0, 3, 3))
  theta <- theta[rep(seq_len(3), c(400, 300, 500)), ]
  y_true <- rowSums(x * theta)
  factor <- c(
    2 * stats::rbinom(400, size = 1, prob = 0.95) - 1,
    2 * stats::rbinom(300, size = 1, prob = 0.95) - 1,
    2 * stats::rbinom(500, size = 1, prob = 0.95) - 1
  )
  y <- factor * y_true + stats::rnorm(n)
  data <- cbind.data.frame(y, x)
  huber_threshold <- 1
  huber_loss <- function(data, theta) {
    residual <- data[, 1] - data[, -1, drop = FALSE] %*% theta
    indicator <- abs(residual) <= huber_threshold
    sum(
      residual^2 / 2 * indicator +
        huber_threshold * (
          abs(residual) - huber_threshold / 2
        ) * (1 - indicator)
    )
  }
  huber_loss_gradient <- function(data, theta) {
    residual <- c(data[nrow(data), 1] - data[nrow(data), -1] %*% theta)
    if (abs(residual) <= huber_threshold) {
      -residual * data[nrow(data), -1]
    } else {
      -huber_threshold * sign(residual) * data[nrow(data), -1]
    }
  }
  huber_loss_hessian <- function(data, theta) {
    residual <- c(data[nrow(data), 1] - data[nrow(data), -1] %*% theta)
    if (abs(residual) <= huber_threshold) {
      outer(data[nrow(data), -1], data[nrow(data), -1])
    } else {
      0.01 * diag(length(theta))
    }
  }
  huber_regression_result <- fastcpd(
    formula = y ~ . - 1,
    data = data,
    beta = (p + 1) * log(n) / 2,
    cost = huber_loss,
    cost_gradient = huber_loss_gradient,
    cost_hessian = huber_loss_hessian
  )
  summary(huber_regression_result)
}
#> 
#> Call:
#> fastcpd(formula = y ~ . - 1, data = data, beta = (p + 1) * log(n)/2, 
#>     cost = huber_loss, cost_gradient = huber_loss_gradient, cost_hessian = huber_loss_hessian)
#> 
#> Change points:
#> 418 726 
#> 
#> Parameters:
#>   segment 1 segment 2 segment 3
#> 1         0         0         0
#> 2         0         0         0
#> 3         0         0         0
#> 4         0         0         0
#> 5         0         0         0
# \donttest{
set.seed(1)
p <- 1
x <- matrix(rnorm(375 * p, 0, 1), ncol = p)
theta <- rbind(rnorm(p, 0, 1), rnorm(p, 2, 1))
y <- c(
  rbinom(200, 1, 1 / (1 + exp(-x[1:200, ] %*% theta[1, , drop = FALSE]))),
  rbinom(175, 1, 1 / (1 + exp(-x[201:375, ] %*% theta[2, , drop = FALSE])))
)
data <- data.frame(y = y, x = x)
result_builtin <- suppressWarnings(fastcpd.binomial(data))
logistic_loss <- function(data, theta) {
  x <- data[, -1, drop = FALSE]
  y <- data[, 1]
  u <- x %*% theta
  nll <- -y * u + log(1 + exp(u))
  nll[u > 10] <- -y[u > 10] * u[u > 10] + u[u > 10]
  sum(nll)
}
logistic_loss_gradient <- function(data, theta) {
  x <- data[nrow(data), -1, drop = FALSE]
  y <- data[nrow(data), 1]
  c(-(y - 1 / (1 + exp(-x %*% theta)))) * x
}
logistic_loss_hessian <- function(data, theta) {
  x <- data[nrow(data), -1]
  prob <- 1 / (1 + exp(-x %*% theta))
  (x %o% x) * c((1 - prob) * prob)
}
result_custom <- fastcpd(
  formula = y ~ . - 1,
  data = data,
  epsilon = 1e-5,
  cost = logistic_loss,
  cost_gradient = logistic_loss_gradient,
  cost_hessian = logistic_loss_hessian
)
result_builtin@cp_set
#> [1] 201
result_custom@cp_set
#> [1] 198
# }
# \donttest{
if (requireNamespace("mvtnorm", quietly = TRUE)) {
  set.seed(1)
  n <- 480
  p_true <- 6
  p <- 50
  x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
  theta_0 <- rbind(
    runif(p_true, -5, -2),
    runif(p_true, -3, 3),
    runif(p_true, 2, 5),
    runif(p_true, -5, 5)
  )
  theta_0 <- cbind(theta_0, matrix(0, ncol = p - p_true, nrow = 4))
  y <- c(
    x[1:80, ] %*% theta_0[1, ] + rnorm(80, 0, 1),
    x[81:200, ] %*% theta_0[2, ] + rnorm(120, 0, 1),
    x[201:320, ] %*% theta_0[3, ] + rnorm(120, 0, 1),
    x[321:n, ] %*% theta_0[4, ] + rnorm(160, 0, 1)
  )
  small_lasso_data <- cbind.data.frame(y, x)
  result_no_vp <- fastcpd.lasso(
    small_lasso_data,
    beta = "BIC",
    cost_adjustment = NULL,
    pruning_coef = 0
  )
  summary(result_no_vp)
  result_20_vp <- fastcpd.lasso(
    small_lasso_data,
    beta = "BIC",
    cost_adjustment = NULL,
    vanilla_percentage = 0.2,
    pruning_coef = 0
  )
  summary(result_20_vp)
}
#> 
#> Call:
#> fastcpd.lasso(data = small_lasso_data, beta = "BIC", cost_adjustment = NULL, 
#>     pruning_coef = 0)
#> 
#> Change points:
#> 79 202 325 
#> 
#> Cost values:
#> 60.77458 1313.455 354.2589 144.8481 
#> 
#> Parameters:
#> 50 x 4 sparse Matrix of class "dgCMatrix"
#>          segment 1  segment 2 segment 3    segment 4
#>  [1,] -2.439372877 -1.2267464  4.199448 -2.871700712
#>  [2,] -2.437178776 -0.8646318  3.691595  2.829794775
#>  [3,] -4.145504088  0.6081470  4.426221  3.496099151
#>  [4,] -4.607393338  .          2.250112 -3.063066325
#>  [5,] -4.574367683  .          4.301097  2.983297359
#>  [6,] -3.808565060 -1.3128104  2.871565  4.067838390
#>  [7,]  .            .          .         .          
#>  [8,]  .            .          .         .          
#>  [9,]  .            .          .         .          
#> [10,]  .            .          .         .          
#> [11,]  0.025537460  .          .        -0.009345011
#> [12,]  .            .          .         .          
#> [13,]  .            .          .         .          
#> [14,]  .            .          .         .          
#> [15,]  .            .          .         .          
#> [16,] -0.018010739  .          .         .          
#> [17,]  .            .          .         .          
#> [18,]  .            .          .         .          
#> [19,]  0.019317827  .          .         .          
#> [20,]  .            .          .         .          
#> [21,]  .            .          .         .          
#> [22,] -0.130706614  .          .         .          
#> [23,]  0.203834430  .          .         .          
#> [24,]  .            .          .         .          
#> [25,] -0.160962429  .          .         .          
#> [26,]  .            .          .         .          
#> [27,]  .            .          .         .          
#> [28,]  0.039164525  .          .         .          
#> [29,]  .            .          .         .          
#> [30,]  .            .          .         .          
#> [31,]  0.004143997  .          .         .          
#> [32,]  .            .          .         .          
#> [33,]  .            .          .         .          
#> [34,]  .            .          .         .          
#> [35,]  .            .          .         .          
#> [36,]  0.193496648  .          .         .          
#> [37,]  .            .          .         0.069864843
#> [38,]  .            .          .         .          
#> [39,]  .            .          .         .          
#> [40,]  .            .          .         .          
#> [41,] -0.066588679  .          .         .          
#> [42,]  .            .          .         .          
#> [43,]  .            .          .         .          
#> [44,]  .            .          .         .          
#> [45,]  .            .          .         .          
#> [46,]  0.033341252  .          .         .          
#> [47,]  .            .          .         .          
#> [48,] -0.019108958  .          .         .          
#> [49,]  .            .          .         .          
#> [50,]  .            .          .         .          
#> 
#> Call:
#> fastcpd.lasso(data = small_lasso_data, beta = "BIC", cost_adjustment = NULL, 
#>     vanilla_percentage = 0.2, pruning_coef = 0)
#> 
#> Change points:
#> 80 202 325 
#> 
#> Cost values:
#> 60.43239 1304.071 354.2589 144.8481 
#> 
#> Parameters:
#> 50 x 4 sparse Matrix of class "dgCMatrix"
#>         segment 1  segment 2 segment 3    segment 4
#>  [1,] -2.45267681 -1.2523758  4.199448 -2.871700712
#>  [2,] -2.47231158 -0.8719669  3.691595  2.829794775
#>  [3,] -4.16314116  0.6099911  4.426221  3.496099151
#>  [4,] -4.58258184  .          2.250112 -3.063066325
#>  [5,] -4.61090633  .          4.301097  2.983297359
#>  [6,] -3.78991394 -1.3094173  2.871565  4.067838390
#>  [7,]  .           .          .         .          
#>  [8,]  .           .          .         .          
#>  [9,]  .           .          .         .          
#> [10,]  .           .          .         .          
#> [11,]  0.04910699  .          .        -0.009345011
#> [12,]  .           .          .         .          
#> [13,]  .           .          .         .          
#> [14,]  .           .          .         .          
#> [15,]  .           .          .         .          
#> [16,] -0.08742899  .          .         .          
#> [17,]  .           .          .         .          
#> [18,]  .           .          .         .          
#> [19,]  .           .          .         .          
#> [20,]  .           .          .         .          
#> [21,]  .           .          .         .          
#> [22,] -0.08053771  .          .         .          
#> [23,]  0.20011359  .          .         .          
#> [24,]  .           .          .         .          
#> [25,] -0.13036946  .          .         .          
#> [26,]  .           .          .         .          
#> [27,]  .           .          .         .          
#> [28,]  0.01682241  .          .         .          
#> [29,]  .           .          .         .          
#> [30,]  .           .          .         .          
#> [31,]  .           .          .         .          
#> [32,]  0.04324226  .          .         .          
#> [33,]  .           .          .         .          
#> [34,]  .           .          .         .          
#> [35,]  .           .          .         .          
#> [36,]  0.25883390  .          .         .          
#> [37,]  .           .          .         0.069864843
#> [38,]  .           .          .         .          
#> [39,]  .           .          .         .          
#> [40,]  .           .          .         .          
#> [41,] -0.10817803  .          .         .          
#> [42,]  .           .          .         .          
#> [43,]  .           .          .         .          
#> [44,]  .           .          .         .          
#> [45,]  .           .          .         .          
#> [46,]  0.03152628  .          .         .          
#> [47,]  .           .          .         .          
#> [48,] -0.02641484  .          .         .          
#> [49,]  .           .          .         .          
#> [50,]  .           .          .         .          
# }
```
