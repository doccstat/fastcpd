#' fastcpd: A package for finding change points in an efficient way
#'
#' The fastcpd package provides a function \code{\link{fastcpd}} to find change
#' points in a data set. The function is based on the paper ``Sequential
#' Gradient Descent and Quasi-Newton's Method for Change-Point Analysis'' by
#' Xianyang Zhang and Trisha Dawn.
#'
#' @section References:
#'   Zhang X, Dawn T (2023). ``Sequential Gradient Descent and Quasi-Newton's
#'   Method for Change-Point Analysis.'' In Ruiz, Francisco, Dy, Jennifer,
#'   van de Meent, Jan-Willem (eds.), _Proceedings of The 26th International
#'   Conference on Artificial Intelligence and Statistics_, volume 206 series
#'   Proceedings of Machine Learning Research, 1129-1143.
#'   <https://proceedings.mlr.press/v206/zhang23b.html>.
#'
#' @docType package
#' @name fastcpd
#' @useDynLib fastcpd, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom methods show
NULL
#> NULL


#' Sequential Gradient Descent and Quasi-Newton's Method for Change-Point
#' Analysis
#'
#' @param formula A formula object specifying the model to be fitted. The
#'   optional response variable should be on the left hand side of the formula
#'   while the covariates should be on the right hand side. The intercept term
#'   should be removed from the formula. The response variable is not
#'   necessary if the data considered is not of regression type. For example,
#'   a mean or variance change model does not necessarily have response
#'   variables. By default an intercept column will be added to the data
#'   similar to the \code{lm} function in \proglang{R}. Thus it is suggested
#'   that user should remove the intercept term from the formula by appending
#'   \code{- 1} to the formula. The default formula is suitable for regression
#'   data sets with one-dimensional response variable and the rest being
#'   covariates without intercept. The naming of variables used in the formula
#'   should be consistent with the column names in the data frame provided in
#'   \code{data}.
#' @param data A data frame containing the data to be segmented where each row
#'   denotes each data point. In one-dimensional response variable regression
#'   settings, the first column is the response variable while the rest are
#'   covariates. The response is not necessary in the case of mean change or
#'   variance change, in which case the formula will need to be adjusted
#'   accordingly.
#' @param beta Initial cost value specified in the algorithm in the paper.
#'   For the proper choice of a value, please refer to the paper. If not
#'   specified, BIC criterion is used to obtain a proper value, i.e.,
#'   \code{beta = (p + 1) * log(nrow(data)) / 2}.
#' @param segment_count Number of segments for initial guess. If not specified,
#'   the initial guess on the number of segments is 10.
#' @param trim Trimming for the boundary change points so that a change point
#'   close to the boundary will not be counted as a change point. This
#'   parameter also specifies the minimum distance between two change points.
#'   If.   several change points have mutual distances smaller than
#'   \code{trim * nrow(data)}, those change points will be merged into one
#'   single change point. The value of this parameter should be between
#'   0 and 1.
#' @param momentum_coef Momentum coefficient to be applied to each update. This
#'   parameter is used when the loss function is bad-shaped so that
#'   maintaining a momentum from previous update is desired. Default value is
#'   0, meaning the algorithm doesn't maintain a momentum by default.
#' @param k Function on number of epochs in SGD. \code{k} should be a function
#'   taking only a parameter \code{x} meaning the current number of data
#'   points considered since last segmentaion. The return value of the
#'   function should be an integer indicating how many epochs should be
#'   performed apart from the default update. By default the function returns
#'   0, meaning no multiple epochs will be used to update the parameters.
#'   Example usage:
#'   ```r
#'     k = function(x) {
#'       if (x < n / segment_count / 4 * 1) 3
#'       else if (x < n / segment_count / 4 * 2) 2
#'       else if (x < n / segment_count / 4 * 3) 1
#'       else 0
#'     }
#'   ```
#'   This function will perform 3 epochs for the first quarter of the data, 2
#'   epochs for the second quarter of the data, 1 epoch for the third quarter
#'   of the data and no multiple epochs for the last quarter of the data.
#'   Experiments show that performing multiple epochs will significantly
#'   affect the performance of the algorithm. This parameter is left for the
#'   users to tune the performance of the algorithm if the result is not
#'   ideal. Details are discussed in the paper.
#' @param family Family of the model. Can be \code{"gaussian"},
#'   \code{"binomial"}, \code{"poisson"}, \code{"lasso"}, \code{"custom"} or
#'   \code{NULL}. For simplicity, user can also omit this parameter,
#'   indicating that they will be using their own cost functions. Omitting the
#'   parameter is the same as specifying the parameter to be \code{"custom"}
#'   or \code{NULL}, in which case, users must specify the cost function, with
#'   optional gradient and corresponding Hessian matrix functions.
#' @param epsilon Epsilon to avoid numerical issues. Only used for the Hessian
#'   computation in Logistic Regression and Poisson Regression.
#' @param min_prob Minimum probability to avoid numerical issues. Only used
#'   for Poisson Regression.
#' @param winsorise_minval Minimum value for the parameter in Poisson Regression
#'   to be winsorised.
#' @param winsorise_maxval Maximum value for the parameter in Poisson Regression
#'   to be winsorised.
#' @param p Number of covariates in the model. If not specified, the number of
#'   covariates will be inferred from the data, i.e.,
#'   \code{p = ncol(data) - 1}.
#' @param cost Cost function to be used. This and the following two parameters
#'   should not be specified at the same time with \code{family}. If not
#'   specified, the default is the negative log-likelihood for the
#'   corresponding family. Custom cost functions can be provided in the
#'   following two formats:
#'
#'   - \code{cost = function(data) \{...\}}
#'   - \code{cost = function(data, theta) \{...\}}
#'
#'   In both methods, users should implement the cost value calculation based
#'   on the data provided, where the data parameter can be considered as a
#'   segment of the original data frame in the form of a matrix. The first
#'   method is used when the cost function has an explicit solution, in which
#'   case the cost function value can be calculated directly from the data.
#'   The second method is used when the cost function does not have an
#'   explicit solution, in which case the cost function value can be
#'   calculated from the data and the estimated parameters. In the case of
#'   only one \code{data} argument is provided, `fastcpd` performs the
#'   vanilla PELT algorithm since no parameter updating is performed.
#' @param cost_gradient Gradient function for the custom cost function.
#'   Example usage:
#'   ```r
#'     cost_gradient = function(data, theta) {
#'       ...
#'       return(gradient)
#'     }
#'   ```
#'   The gradient function should take two parameters, the first one being a
#'   segment of the data in the format of a matrix, the second one being the
#'   estimated parameters. The gradient function should return the gradient of
#'   the cost function with respect to the data and parameters.
#' @param cost_hessian Hessian function for the custom cost function. Similar to
#'   the gradient function, the Hessian function should take two parameters,
#'   the first one being a segment of the data in the format of a matrix, the
#'   second one being the estimated parameters. The Hessian function should
#'   return the Hessian matrix of the cost function with respect to the data
#'   and parameters.
#' @param cp_only If \code{TRUE}, only the change points are returned.
#'   Otherwise, the cost function values together with the estimated
#'   parameters for each segment are also returned. By default the value is
#'   set to be \code{FALSE} so that `plot` can be used to visualize the
#'   results for a built-in model. If \code{family} is not provided or
#'   specified as \code{NULL} or \code{"custom"}, \code{cp_only} is set to be
#'   \code{TRUE} by default. \code{cp_only} has some performance impact on the
#'   algorithm, since the cost values and estimated parameters for each
#'   segment need to be calculated and stored. If the users are only
#'   interested in the change points, setting \code{cp_only} to be \code{TRUE}
#'   will help with the computational cost.
#' @param vanilla_percentage How many of the data should be processed through
#'   vanilla PELT. Range should be between 0 and 1. The `fastcpd`
#'   algorithm is based on gradient descent and thus a starting estimate can
#'   be crucial. At the beginning of the algorithm, vanilla PELT can be
#'   performed to obtain a relatively accurate estimate of the parameters
#'   despite the small amount of the data being used. If set to be 0, all data
#'   will be processed through sequential gradient descnet. If set to be 1,
#'   all data will be processed through vaniall PELT. If the cost function
#'   have an explicit solution, i.e. does not depend on coefficients like the
#'   mean change case, this parameter will be set to be 1. If the value is set
#'   to be between 0 and 1, the first \code{vanilla_percentage * nrow(data)}
#'   data points will be processed through vanilla PELT and the rest will be
#'   processed through sequential gradient descent.
#'
#' @return A class \code{fastcpd} object.
#' @export
#' @md
#' @examples
#' \donttest{
#' ### linear regression
#' library(fastcpd)
#' set.seed(1)
#' p <- 3
#' x <- mvtnorm::rmvnorm(300, rep(0, p), diag(p))
#' theta_0 <- rbind(c(1, 1.2, -1), c(-1, 0, 0.5), c(0.5, -0.3, 0.2))
#' y <- c(
#'   x[1:100, ] %*% theta_0[1, ] + rnorm(100, 0, 1),
#'   x[101:200, ] %*% theta_0[2, ] + rnorm(100, 0, 1),
#'   x[201:300, ] %*% theta_0[3, ] + rnorm(100, 0, 1)
#' )
#' result <- fastcpd(
#'   formula = y ~ . - 1,
#'   data = data.frame(y = y, x = x),
#'   family = "gaussian"
#' )
#' plot(result)
#' summary(result)
#'
#' ### linear regression with one-dimensional covariate
#' library(fastcpd)
#' set.seed(1)
#' p <- 1
#' x <- mvtnorm::rmvnorm(300, rep(0, p), diag(p))
#' theta_0 <- matrix(c(1, -1, 0.5))
#' y <- c(
#'   x[1:100, ] * theta_0[1, ] + rnorm(100, 0, 1),
#'   x[101:200, ] * theta_0[2, ] + rnorm(100, 0, 1),
#'   x[201:300, ] * theta_0[3, ] + rnorm(100, 0, 1)
#' )
#' result <- fastcpd(
#'   formula = y ~ . - 1,
#'   data = data.frame(y = y, x = x),
#'   family = "gaussian"
#' )
#' plot(result)
#' summary(result)
#'
#' ### logistic regression
#' library(fastcpd)
#' set.seed(1)
#' x <- matrix(rnorm(1500, 0, 1), ncol = 5)
#' theta <- rbind(rnorm(5, 0, 1), rnorm(5, 2, 1))
#' y <- c(
#'   rbinom(125, 1, 1 / (1 + exp(-x[1:125, ] %*% theta[1, ]))),
#'   rbinom(175, 1, 1 / (1 + exp(-x[126:300, ] %*% theta[2, ])))
#' )
#' result <- suppressWarnings(fastcpd(
#'   formula = y ~ . - 1,
#'   data = data.frame(y = y, x = x),
#'   family = "binomial"
#' ))
#' summary(result)
#'
#' ### poisson regression
#' library(fastcpd)
#' set.seed(1)
#' p <- 3
#' x <- mvtnorm::rmvnorm(1500, rep(0, p), diag(p))
#' delta <- rnorm(p)
#' theta_0 <- c(1, 1.2, -1)
#' y <- c(
#'   rpois(300, exp(x[1:300, ] %*% theta_0)),
#'   rpois(400, exp(x[301:700, ] %*% (theta_0 + delta))),
#'   rpois(300, exp(x[701:1000, ] %*% theta_0)),
#'   rpois(100, exp(x[1001:1100, ] %*% (theta_0 - delta))),
#'   rpois(200, exp(x[1101:1300, ] %*% theta_0)),
#'   rpois(200, exp(x[1301:1500, ] %*% (theta_0 + delta)))
#' )
#' result <- fastcpd(
#'   formula = y ~ . - 1,
#'   data = data.frame(y = y, x = x),
#'   beta = (p + 1) * log(1500) / 2,
#'   k = function(x) 0,
#'   family = "poisson",
#'   epsilon = 1e-5
#' )
#' summary(result)
#' result_two_epochs <- fastcpd(
#'   formula = y ~ . - 1,
#'   data = data.frame(y = y, x = x),
#'   beta = (p + 1) * log(1500) / 2,
#'   k = function(x) 1,
#'   family = "poisson",
#'   epsilon = 1e-4
#' )
#' summary(result_two_epochs)
#'
#' ### penalized linear regression
#' library(fastcpd)
#' set.seed(1)
#' n <- 1500
#' p_true <- 6
#' p <- 50
#' x <- mvtnorm::rmvnorm(1500, rep(0, p), diag(p))
#' theta_0 <- rbind(
#'   runif(p_true, -5, -2),
#'   runif(p_true, -3, 3),
#'   runif(p_true, 2, 5),
#'   runif(p_true, -5, 5)
#' )
#' theta_0 <- cbind(theta_0, matrix(0, ncol = p - p_true, nrow = 4))
#' y <- c(
#'   x[1:300, ] %*% theta_0[1, ] + rnorm(300, 0, 1),
#'   x[301:700, ] %*% theta_0[2, ] + rnorm(400, 0, 1),
#'   x[701:1000, ] %*% theta_0[3, ] + rnorm(300, 0, 1),
#'   x[1001:1500, ] %*% theta_0[4, ] + rnorm(500, 0, 1)
#' )
#' result <- fastcpd(
#'   formula = y ~ . - 1,
#'   data = data.frame(y = y, x = x),
#'   family = "lasso"
#' )
#' plot(result)
#' summary(result)
#'
#' ### custom logistic regression
#' library(fastcpd)
#' set.seed(1)
#' p <- 5
#' x <- matrix(rnorm(375 * p, 0, 1), ncol = p)
#' theta <- rbind(rnorm(p, 0, 1), rnorm(p, 2, 1))
#' y <- c(
#'   rbinom(200, 1, 1 / (1 + exp(-x[1:200, ] %*% theta[1, ]))),
#'   rbinom(175, 1, 1 / (1 + exp(-x[201:375, ] %*% theta[2, ])))
#' )
#' data <- data.frame(y = y, x = x)
#' result_builtin <- suppressWarnings(fastcpd(
#'   formula = y ~ . - 1,
#'   data = data,
#'   family = "binomial"
#' ))
#' logistic_loss <- function(data, theta) {
#'   x <- data[, -1]
#'   y <- data[, 1]
#'   u <- x %*% theta
#'   nll <- -y * u + log(1 + exp(u))
#'   nll[u > 10] <- -y[u > 10] * u[u > 10] + u[u > 10]
#'   sum(nll)
#' }
#' logistic_loss_gradient <- function(data, theta) {
#'   x <- data[nrow(data), -1]
#'   y <- data[nrow(data), 1]
#'   c(-(y - 1 / (1 + exp(-x %*% theta)))) * x
#' }
#' logistic_loss_hessian <- function(data, theta) {
#'   x <- data[nrow(data), -1]
#'   prob <- 1 / (1 + exp(-x %*% theta))
#'   (x %o% x) * c((1 - prob) * prob)
#' }
#' result_custom <- fastcpd(
#'   formula = y ~ . - 1,
#'   data = data,
#'   epsilon = 1e-5,
#'   cost = logistic_loss,
#'   cost_gradient = logistic_loss_gradient,
#'   cost_hessian = logistic_loss_hessian
#' )
#' cat(
#'   "Change points detected by built-in logistic regression model: ",
#'   result_builtin@cp_set, "\n",
#'   "Change points detected by custom logistic regression model: ",
#'   result_custom@cp_set, "\n",
#'   sep = ""
#' )
#' result_custom_two_epochs <- fastcpd(
#'   formula = y ~ . - 1,
#'   data = data,
#'   k = function(x) 1,
#'   epsilon = 1e-5,
#'   cost = logistic_loss,
#'   cost_gradient = logistic_loss_gradient,
#'   cost_hessian = logistic_loss_hessian
#' )
#' summary(result_custom_two_epochs)
#'
#' ### custom cost function mean change
#' library(fastcpd)
#' set.seed(1)
#' p <- 1
#' data <- rbind(
#'   mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(100, p)),
#'   mvtnorm::rmvnorm(400, mean = rep(50, p), sigma = diag(100, p)),
#'   mvtnorm::rmvnorm(300, mean = rep(2, p), sigma = diag(100, p))
#' )
#' segment_count_guess <- 10
#' block_size <- max(floor(sqrt(nrow(data)) / (segment_count_guess + 1)), 2)
#' block_count <- floor(nrow(data) / block_size)
#' data_all_vars <- rep(0, block_count)
#' for (block_index in seq_len(block_count)) {
#'   block_start <- (block_index - 1) * block_size + 1
#'   block_end <- if (block_index < block_count) {
#'     block_index * block_size
#'   } else {
#'     nrow(data)
#'   }
#'   data_all_vars[block_index] <- var(data[block_start:block_end, ])
#' }
#' data_all_var <- mean(data_all_vars)
#' mean_loss <- function(data) {
#'   n <- nrow(data)
#'   n / 2 * (
#'     log(data_all_var) + log(2 * pi) +
#'       sum((data - colMeans(data))^2 / data_all_var) / n
#'   )
#' }
#' mean_loss_result <- fastcpd(
#'   formula = ~ . - 1,
#'   data = data.frame(data),
#'   beta = (p + 1) * log(nrow(data)) / 2,
#'   p = p,
#'   cost = mean_loss
#' )
#' summary(mean_loss_result)
#'
#' ### custom cost function multivariate mean change
#' library(fastcpd)
#' set.seed(1)
#' p <- 3
#' data <- rbind(
#'   mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(100, p)),
#'   mvtnorm::rmvnorm(400, mean = rep(50, p), sigma = diag(100, p)),
#'   mvtnorm::rmvnorm(300, mean = rep(2, p), sigma = diag(100, p))
#' )
#' segment_count_guess <- 5
#' block_size <- max(floor(sqrt(nrow(data)) / (segment_count_guess + 1)), 2)
#' block_count <- floor(nrow(data) / block_size)
#' data_all_covs <- array(NA, dim = c(block_count, p, p))
#' for (block_index in seq_len(block_count)) {
#'   block_start <- (block_index - 1) * block_size + 1
#'   block_end <- if (block_index < block_count) {
#'     block_index * block_size
#'   } else {
#'     nrow(data)
#'   }
#'   data_all_covs[block_index, , ] <- cov(data[block_start:block_end, ])
#' }
#' data_all_cov <- colMeans(data_all_covs)
#' mean_loss <- function(data) {
#'   n <- nrow(data)
#'   demeaned_data <- sweep(data, 2, colMeans(data))
#'   n / 2 * (
#'     log(det(data_all_cov)) + p * log(2 * pi) +
#'       sum(diag(solve(data_all_cov, crossprod(demeaned_data)))) / n
#'   )
#' }
#' mean_loss_result <- fastcpd(
#'   formula = ~ . - 1,
#'   data = data.frame(data),
#'   beta = (p + 1) * log(nrow(data)) / 2,
#'   p = p,
#'   cost = mean_loss
#' )
#' summary(mean_loss_result)
#'
#' ### custom cost function variance change
#' library(fastcpd)
#' set.seed(1)
#' p <- 1
#' data <- rbind.data.frame(
#'   mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
#'   mvtnorm::rmvnorm(400, mean = rep(0, p), sigma = diag(50, p)),
#'   mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(2, p))
#' )
#' data_all_mean <- colMeans(data)
#' var_loss <- function(data) {
#'   n <- nrow(data)
#'   data_cov <- crossprod(sweep(data, 2, data_all_mean)) / (n - 1)
#'   n / 2 * (log(data_cov) + log(2 * pi) + (n - 1) / n)
#' }
#' var_loss_result <- fastcpd(
#'   formula = ~ . - 1,
#'   data = data,
#'   beta = (p + 1) * log(nrow(data)) / 2,
#'   p = p,
#'   cost = var_loss
#' )
#' summary(var_loss_result)
#'
#' ### custom cost function multivariate variance change
#' library(fastcpd)
#' set.seed(1)
#' p <- 3
#' data <- rbind.data.frame(
#'   mvtnorm::rmvnorm(
#'     300, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))
#'   ),
#'   mvtnorm::rmvnorm(
#'     400, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))
#'   ),
#'   mvtnorm::rmvnorm(
#'     300, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))
#'   )
#' )
#' data_all_mean <- colMeans(data)
#' var_loss <- function(data) {
#'   n <- nrow(data)
#'   p <- ncol(data)
#'   if (n < p) {
#'     data_cov <- diag(p)
#'   } else {
#'     data_cov <- crossprod(sweep(data, 2, data_all_mean)) / (n - 1)
#'   }
#'   n / 2 * (log(det(data_cov)) + p * log(2 * pi) + p * (n - 1) / n)
#' }
#' var_loss_result <- fastcpd(
#'   formula = ~ . - 1,
#'   data = data,
#'   beta = (p^2 + 1) * log(nrow(data)) / 2,
#'   trim = 0.1,
#'   p = p^2,
#'   cost = var_loss
#' )
#' summary(var_loss_result)
#'
#' ### custom cost function mean or variance change
#' library(fastcpd)
#' set.seed(1)
#' p <- 1
#' data <- rbind.data.frame(
#'   mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
#'   mvtnorm::rmvnorm(400, mean = rep(10, p), sigma = diag(1, p)),
#'   mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(50, p)),
#'   mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
#'   mvtnorm::rmvnorm(400, mean = rep(10, p), sigma = diag(1, p)),
#'   mvtnorm::rmvnorm(300, mean = rep(10, p), sigma = diag(50, p))
#' )
#' meanvar_loss <- function(data) {
#'   n <- nrow(data)
#'   data_cov <- 1
#'   if (n > 1) {
#'     data_cov <- var(data)
#'   }
#'   n / 2 * (log(data_cov) + log(2 * pi) + (n - 1) / n)
#' }
#' meanvar_loss_result <- fastcpd(
#'   formula = ~ . - 1,
#'   data = data,
#'   beta = (p^2 + p + 1) * log(nrow(data)) / 2,
#'   p = p^2 + p,
#'   cost = meanvar_loss
#' )
#' summary(meanvar_loss_result)
#'
#' ### custom cost function multivariate mean or variance change
#' library(fastcpd)
#' set.seed(1)
#' p <- 3
#' data <- rbind.data.frame(
#'   mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
#'   mvtnorm::rmvnorm(400, mean = rep(10, p), sigma = diag(1, p)),
#'   mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(50, p)),
#'   mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
#'   mvtnorm::rmvnorm(400, mean = rep(10, p), sigma = diag(1, p)),
#'   mvtnorm::rmvnorm(300, mean = rep(10, p), sigma = diag(50, p))
#' )
#' meanvar_loss <- function(data) {
#'   n <- nrow(data)
#'   p <- ncol(data)
#'   if (n <= p) {
#'     data_cov <- diag(p)
#'   } else {
#'     data_cov <- cov(data)
#'   }
#'   n / 2 * (log(det(data_cov)) + p * log(2 * pi) + p * (n - 1) / n)
#' }
#' meanvar_loss_result <- fastcpd(
#'   formula = ~ . - 1,
#'   data = data,
#'   beta = (p^2 + p + 1) * log(nrow(data)) / 2,
#'   trim = 0.01,
#'   p = p^2 + p,
#'   cost = meanvar_loss
#' )
#' summary(meanvar_loss_result)
#'
#' ### custom cost function huber regression
#' library(fastcpd)
#' set.seed(1)
#' n <- 400 + 300 + 500
#' p <- 5
#' x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
#' theta <- rbind(
#'   mvtnorm::rmvnorm(1, mean = rep(0, p - 3), sigma = diag(p - 3)),
#'   mvtnorm::rmvnorm(1, mean = rep(5, p - 3), sigma = diag(p - 3)),
#'   mvtnorm::rmvnorm(1, mean = rep(9, p - 3), sigma = diag(p - 3))
#' )
#' theta <- cbind(theta, matrix(0, 3, 3))
#' theta <- theta[rep(seq_len(3), c(400, 300, 500)), ]
#' y_true <- rowSums(x * theta)
#' factor <- c(
#'   2 * stats::rbinom(400, size = 1, prob = 0.95) - 1,
#'   2 * stats::rbinom(300, size = 1, prob = 0.95) - 1,
#'   2 * stats::rbinom(500, size = 1, prob = 0.95) - 1
#' )
#' y <- factor * y_true + stats::rnorm(n)
#' data <- cbind.data.frame(y, x)
#' huber_threshold <- 1
#' huber_loss <- function(data, theta) {
#'   residual <- data[, 1] - data[, -1, drop = FALSE] %*% theta
#'   indicator <- abs(residual) <= huber_threshold
#'   sum(
#'     residual^2 / 2 * indicator +
#'       huber_threshold * (
#'         abs(residual) - huber_threshold / 2
#'       ) * (1 - indicator)
#'   )
#' }
#' huber_loss_gradient <- function(data, theta) {
#'   residual <- c(data[nrow(data), 1] - data[nrow(data), -1] %*% theta)
#'   if (abs(residual) <= huber_threshold) {
#'     -residual * data[nrow(data), -1]
#'   } else {
#'     -huber_threshold * sign(residual) * data[nrow(data), -1]
#'   }
#' }
#' huber_loss_hessian <- function(data, theta) {
#'   residual <- c(data[nrow(data), 1] - data[nrow(data), -1] %*% theta)
#'   if (abs(residual) <= huber_threshold) {
#'     outer(data[nrow(data), -1], data[nrow(data), -1])
#'   } else {
#'     0.01 * diag(length(theta))
#'   }
#' }
#' huber_regression_result <- fastcpd(
#'   formula = y ~ . - 1,
#'   data = data,
#'   beta = (p + 1) * log(n) / 2,
#'   cost = huber_loss,
#'   cost_gradient = huber_loss_gradient,
#'   cost_hessian = huber_loss_hessian
#' )
#' summary(huber_regression_result)
#' }
fastcpd <- function(
    formula = y ~ . - 1,
    data,
    beta = NULL,
    segment_count = 10,
    trim = 0.025,
    momentum_coef = 0,
    k = function(x) 0,
    family = NULL,
    epsilon = 1e-10,
    min_prob = 10^10,
    winsorise_minval = -20,
    winsorise_maxval = 20,
    p = NULL,
    cost = negative_log_likelihood,
    cost_gradient = cost_update_gradient,
    cost_hessian = cost_update_hessian,
    cp_only = FALSE,
    vanilla_percentage = 0) {
  # The following code is adapted from the `lm` function from base R.
  match_formula <- match.call(expand.dots = FALSE)
  matched_formula <- match(c("formula", "data"), names(match_formula), 0L)
  match_formula <- match_formula[c(1L, matched_formula)]
  match_formula$drop.unused.levels <- TRUE
  match_formula[[1L]] <- quote(stats::model.frame)
  match_formula <- eval(match_formula, parent.frame())
  mt <- attr(match_formula, "terms")
  y <- stats::model.response(match_formula, "numeric")
  x <- stats::model.matrix(mt, match_formula)
  data <- cbind(y, x)

  if (is.null(family) || length(formals(cost)) == 1) {
    family <- "custom"
  }

  if (family == "custom") {
    cp_only <- TRUE
  }

  if (is.null(p)) {
    p <- ncol(data) - 1
  }

  if (is.null(beta)) {
    beta <- (p + 1) * log(nrow(data)) / 2
  }

  if (length(formals(cost)) == 1) {
    vanilla_percentage <- 1
  }

  if (family == "lasso") {
    # TODO(doccstat): Due to the excessive calls to `glmnet` between R and C++,
    # it is better to use the R implementation of `fastcpd` for lasso.
    result <- fastcpd_impl_r(
      data, beta, segment_count, trim, momentum_coef, k, family, epsilon,
      min_prob, winsorise_minval, winsorise_maxval, p, cost, cost_gradient,
      cost_hessian, cp_only, vanilla_percentage
    )
  } else {
    result <- fastcpd_impl(
      data, beta, segment_count, trim, momentum_coef, k, family, epsilon,
      min_prob, winsorise_minval, winsorise_maxval, p, cost, cost_gradient,
      cost_hessian, cp_only, vanilla_percentage
    )
  }

  result$thetas <- data.frame(result$thetas)
  if (ncol(result$thetas) > 0) {
    names(result$thetas) <- paste0("segment ", seq_len(ncol(result$thetas)))
  }

  if (is.null(result$cost_values)) {
    result$cost_values <- numeric(0)
  }

  if (is.null(result$residual)) {
    result$residual <- numeric(0)
  }

  methods::new(
    Class = "fastcpd",
    call = match.call(),
    data = data.frame(data),
    family = family,
    cp_set = c(result$cp_set),
    cost_values = c(result$cost_values),
    residuals = c(result$residual),
    thetas = result$thetas,
    cp_only = cp_only
  )
}

fastcpd_impl_r <- function(
    data, beta, segment_count, trim, momentum_coef, k, family, epsilon,
    min_prob, winsorise_minval, winsorise_maxval, p, cost, cost_gradient,
    cost_hessian, cp_only, vanilla_percentage) {
  # Set up the initial values.
  n <- nrow(data)
  lambda <- 0

  # After t = 1, the r_t_set R_t contains 0 and 1.
  r_t_set <- c(0, 1)
  # C(0)=NULL, C(1)={0}
  cp_set <- append(list(NULL), rep(list(0), n))
  # Objective function: F(0) = -beta
  f_t <- c(-beta, rep(0, n))

  fastcpd_parameters <- init_fastcpd_parameters(
    data, p, family, segment_count, cost, winsorise_minval, winsorise_maxval,
    epsilon, vanilla_percentage, beta
  )

  for (t in 2:n) {
    r_t_count <- length(r_t_set)

    # number of cost values is the same as number of elemnts in R_t
    cval <- rep(0, r_t_count)

    # for tau in R_t\{t-1}
    for (i in 1:(r_t_count - 1)) {
      tau <- r_t_set[i]
      if (family == "lasso") {
        # Mean of `err_sd` only works if error sd is unchanged.
        lambda <- mean(fastcpd_parameters$err_sd) * sqrt(2 * log(p) / (t - tau))
      }

      data_segment <- data[(tau + 1):t, , drop = FALSE]
      if (vanilla_percentage == 1 || t <= vanilla_percentage * n) {
        cost_optim_result <- cost_optim(family, p, data_segment, cost, 0, TRUE)
        cval[i] <- cost_optim_result$value
        if (vanilla_percentage < 1 && t == vanilla_percentage * n) {
          fastcpd_parameters$theta_hat[, i] <- cost_optim_result$par
          fastcpd_parameters$theta_sum[, i] <-
            fastcpd_parameters$theta_sum[, i] + cost_optim_result$par
        }
      } else {
        fastcpd_parameters <- update_fastcpd_parameters(
          fastcpd_parameters, data, t, i, k, tau, lambda, family,
          cost_gradient, cost_hessian, r_t_set, p,
          momentum_coef, min_prob, winsorise_minval, winsorise_maxval, epsilon
        )

        theta <- fastcpd_parameters$theta_sum[, i] / (t - tau)
        if (family == "poisson" && t - tau >= p) {
          theta <- DescTools::Winsorize(
            x = theta, minval = winsorise_minval, maxval = winsorise_maxval
          )
        }

        if (
          (family %in% c("gaussian", "binomial", "poisson") && t - tau >= p) ||
            (family == "lasso" && t - tau >= 3)
        ) {
          cval[i] <- cost(data_segment, theta, family, lambda)$value
        } else if (family == "custom") {
          # if (warm_start && t - tau >= 50) {
          #   cost_result <- cost(data_segment, start = start[, tau + 1])
          #   start[, tau + 1] <- cost_result$par
          #   cval[i] <- cost_result$value
          # } else {
          cval[i] <- cost(data_segment, theta)
          # }
        }
      }
    }
    fastcpd_parameters <- append_fastcpd_parameters(
      fastcpd_parameters, vanilla_percentage, data, t, family,
      winsorise_minval, winsorise_maxval, p, epsilon
    )

    # Step 3
    cval[r_t_count] <- 0
    # `beta` adjustment seems to work but there might be better choices.
    obj <- cval + f_t[r_t_set + 1] + beta
    min_obj <- min(obj)
    tau_star <- r_t_set[which(obj == min_obj)[1]]

    # Step 4
    cp_set[[t + 1]] <- c(cp_set[[tau_star + 1]], tau_star)

    # Step 5
    pruned_left <- (cval + f_t[r_t_set + 1]) <= min_obj
    r_t_set <- c(r_t_set[pruned_left], t)

    if (vanilla_percentage != 1) {
      fastcpd_parameters$theta_hat <-
        fastcpd_parameters$theta_hat[, pruned_left, drop = FALSE]
      fastcpd_parameters$theta_sum <-
        fastcpd_parameters$theta_sum[, pruned_left, drop = FALSE]
      fastcpd_parameters$hessian <-
        fastcpd_parameters$hessian[, , pruned_left, drop = FALSE]
    }

    # Objective function F(t).
    f_t[t + 1] <- min_obj
  }

  # Remove change-points close to the boundaries
  cp_set <- cp_set[[n + 1]]
  cp_set <- cp_set[(cp_set >= trim * n) & (cp_set <= (1 - trim) * n)]
  cp_set <- sort(unique(c(0, cp_set)))

  cp_set_too_close <- which((diff(cp_set) < trim * n) == TRUE)
  if (length(cp_set_too_close) > 0) {
    cp_set <- floor(
      (cp_set[-(cp_set_too_close + 1)] + cp_set[-cp_set_too_close]) / 2
    )
  }
  cp_set <- cp_set[cp_set > 0]

  residual <- numeric(0)
  if (cp_only) {
    cost_values <- numeric(0)
    thetas <- matrix(NA, nrow = 0, ncol = 0)
  } else {
    cp_loc <- unique(c(0, cp_set, n))
    cost_values <- rep(0, length(cp_loc) - 1)
    thetas <- matrix(NA, nrow = p, ncol = length(cp_loc) - 1)
    for (i in 1:(length(cp_loc) - 1)) {
      data_segment <- data[(cp_loc[i] + 1):cp_loc[i + 1], , drop = FALSE]
      cost_result <- cost_optim(family, p, data_segment, cost, lambda, FALSE)
      residual <- c(residual, cost_result$residuals)
      cost_values[i] <- cost_result$value
      thetas[, i] <- cost_result$par
    }
  }
  thetas <- data.frame(thetas)
  if (ncol(thetas) > 0) {
    names(thetas) <- paste0("segment ", seq_len(ncol(thetas)))
  }
  list(
    cp_set = cp_set,
    cost_values = cost_values,
    residual = residual,
    thetas = thetas
  )
}
