#' @title Find change points efficiently
#' @param formula A formula object specifying the model to be fitted. The
#' (optional) response variable should be on the LHS of the formula, while the
#' covariates should be on the RHS. The naming of variables used in the formula
#' should be consistent with the column names in the data frame provided in
#' \code{data}. The intercept term should be removed from the formula.
#' The response variable is not needed for mean/variance change models and time
#' series models. By default, an intercept column will be added to the data,
#' similar to the [lm()] function.
#' Thus, it is suggested that users should remove the intercept term by
#' appending \code{- 1} to the formula. Note that the [fastcpd.family] functions
#' do not require a formula input.
#' @param data A data frame of dimension \eqn{T \times d}{T * d} containing the
#' data to be segmented (where each row denotes a data point
#' \eqn{z_t \in \mathbb{R}^d}{z_t in R^d} for \eqn{t = 1, \ldots, T}) is
#' required in the main function, while a matrix or a vector input is also
#' accepted in the [fastcpd.family] functions.
#' @param beta Penalty criterion for the number of change points. This parameter
#' takes a string value of \code{"BIC"}, \code{"MBIC"}, \code{"MDL"} or a
#' numeric value.
#' If a numeric value is provided, the value will be used as the penalty.
#' By default, the mBIC criterion is used, where
#' \eqn{\beta = (p + 2) \log(T) / 2}{\beta = (p + 2) log(T) / 2}.
#' This parameter usage should be paired with \code{cost_adjustment} described
#' below. Discussions about the penalty criterion can be found in the
#' references.
#' @param cost_adjustment Cost adjustment criterion.
#' It can be \code{"BIC"}, \code{"MBIC"}, \code{"MDL"} or \code{NULL}.
#' By default, the cost adjustment criterion is set to be \code{"MBIC"}.
#' The \code{"MBIC"} and \code{"MDL"} criteria modify the cost function by
#' adding a negative adjustment term to the cost function.
#' \code{"BIC"} or \code{NULL} does not modify the cost function.
#' Details can in found in the references.
#' @param family Family class of the change point model. It can be
#' \code{"mean"} for mean change,
#' \code{"variance"} for variance change,
#' \code{"meanvariance"} for mean and/or variance change,
#' \code{"lm"} for linear regression,
#' \code{"binomial"} for logistic regression,
#' \code{"poisson"} for Poisson regression,
#' \code{"lasso"} for penalized linear regression,
#' \code{"ar"} for AR(\eqn{p}) models,
#' \code{"arma"} for ARMA(\eqn{p}, \eqn{q}) models,
#' \code{"arima"} for ARIMA(\eqn{p}, \eqn{d}, \eqn{q}) models,
#' \code{"garch"} for GARCH(\eqn{p}, \eqn{q}) models,
#' \code{"var"} for VAR(\eqn{p}) models and
#' \code{"custom"} for user-specified custom models.
#' Omitting this parameter is the same as specifying the parameter to be
#' \code{"custom"} or \code{NULL}, in which case, users must specify the
#' custom cost function.
#' @param cost Cost function to be used. \code{cost}, \code{cost_gradient}, and
#' \code{cost_hessian} should not be specified at the same time with
#' \code{family} as built-in families have cost functions implemented in C++
#' to provide better performance. If not specified, the default is the negative
#' log-likelihood for the corresponding family. Custom cost functions can be
#' provided in the following two formats:
#' \itemize{
#' \item \code{cost = function(data) \{...\}}
#' \item \code{cost = function(data, theta) \{...\}}
#' }
#' Users can specify a loss function using the second format that will be used
#' to calculate the cost value. In both formats, the input data is a subset of
#' the original data frame in the form of a matrix
#' (a matrix with a single column in the case of a univariate data set).
#' In the first format, the specified cost function directly calculates the cost
#' value. [fastcpd()] performs the vanilla PELT algorithm, and
#' \code{cost_gradient} and \code{cost_hessian} should not be provided since no
#' parameter updating is necessary for vanilla PELT.
#' In the second format, the loss function
#' \eqn{\sum_{i = s}^t l(z_i, \theta)}{sum_{i = s}^t l(z_i, \theta)} is
#' provided, which has to be optimized over the parameter \eqn{\theta} to
#' obtain the cost value. A detailed discussion about the custom cost function
#' usage can be found in the references.
#' @param cost_gradient Gradient of the custom cost function. Example usage:
#' ```r
#' cost_gradient = function(data, theta) {
#'   ...
#'   return(gradient)
#' }
#' ```
#' The gradient function takes two inputs, the first being a matrix representing
#' a segment of the data, similar to the format used in the \code{cost}
#' function, and the second being the parameter that needs to be optimized.
#' The gradient function returns the value of the gradient of the loss function,
#' i.e.,
#' \eqn{\sum_{i = s}^t \nabla l(z_i, \theta)}{sum_{i = s}^t l'(z_i, \theta)}.
#' @param cost_hessian Hessian of the custom loss function. The Hessian function
#' takes two inputs, the first being a matrix representing a segment of the
#' data, similar to the format used in the \code{cost} function, and the second
#' being the parameter that needs to be optimized. The gradient function returns
#' the Hessian of the loss function, i.e.,
#' \eqn{\sum_{i = s}^t \nabla^2 l(z_i, \theta)}{sum_{i = s}^t l''(z_i, \theta)}.
#' @param line_search If a vector of numeric values is provided, a line search
#' will be performed to find the optimal step size for each update. Detailed
#' usage of \code{line_search} can be found in the references.
#' @param lower Lower bound for the parameters. Used to specify the domain of
#' the parameters after each gradient descent step. If not specified, the lower
#' bound is set to be \code{-Inf} for all parameters. \code{lower} is especially
#' useful when the estimated parameters take only positive values, such as the
#' noise variance.
#' @param upper Upper bound for the parameters. Used to specify the domain of
#' the parameters after each gradient descent step. If not specified, the upper
#' bound is set to be \code{Inf} for all parameters.
#' @param pruning_coef Pruning coefficient $c_0$ used in the pruning step of the
#' PELT algorithm with the default value 0. If \code{cost_adjustment} is
#' specified as \code{"MBIC"}, an adjustment term \eqn{p\log(2)}{p * log(2)}
#' will be added to the pruning coefficient. If \code{cost_adjustment} is
#' specified as \code{"MDL"}, an adjustment term \eqn{p\log_2(2)}{p * log2(2)}
#' will be added to the pruning coefficient. Detailed discussion about the
#' pruning coefficient can be found in the references.
#' @param segment_count An initial guess of the number of segments. If not
#' specified, the initial guess of the number of segments is 10. The initial
#' guess affects the initial estimates of the parameters in SeGD.
#' @param trim Trimming for the boundary change points so that a change point
#' close to the boundary will not be counted as a change point. This
#' parameter also specifies the minimum distance between two change points.
#' If several change points have mutual distances smaller than
#' \code{trim * nrow(data)}, those change points will be merged into one
#' single change point. The value of this parameter should be between
#' 0 and 1.
#' @param momentum_coef Momentum coefficient to be applied to each update. This
#' parameter is used when the loss function is bad-shaped so that
#' maintaining a momentum from previous update is desired. Default value is
#' 0, meaning the algorithm doesn't maintain a momentum by default.
#' @param multiple_epochs A function can be specified such that an adaptive
#' number of multiple epochs can be utilized to improve the algorithm's
#' performance. \code{multiple_epochs} is a function of the length of the data
#' segment. The function returns an integer indicating how many epochs should be
#' performed apart from the default update. By default, the function returns
#' zero, meaning no multiple epochs will be used to update the parameters.
#' Example usage:
#' ```r
#' multiple_epochs = function(segment_length) {
#'   if (segment_length < 100) 1
#'   else 0
#' }
#' ```
#' This function will let SeGD perform parameter updates with an additional
#' epoch for each segment with a length less than 100 and no additional epoch
#' for segments with lengths greater or equal to 100.
#' @param epsilon Epsilon to avoid numerical issues. Only used for the Hessian
#' computation in Logistic Regression and Poisson Regression.
#' @param order Order of the AR(\eqn{p}), VAR(\eqn{p}) or
#' ARIMA(\eqn{p}, \eqn{d}, \eqn{q}) model.
#' @param p Number of covariates in the model. If not specified, the number of
#' covariates will be inferred from the data, i.e.,
#' \code{p = ncol(data) - 1}. This parameter is superseded by `order` in the
#' case of time series models: "ar", "var", "arima".
#' @param cp_only If \code{TRUE}, only the change points are returned.
#' Otherwise, the cost function values together with the estimated
#' parameters for each segment are also returned. By default the value is
#' set to be \code{FALSE} so that `plot` can be used to visualize the
#' results for a built-in model. \code{cp_only} has some performance impact
#' on the algorithm, since the cost values and estimated parameters for each
#' segment need to be calculated and stored. If the users are only
#' interested in the change points, setting \code{cp_only} to be \code{TRUE}
#' will help with the computational cost.
#' @param vanilla_percentage The parameter \eqn{v} is between zero and one.
#' For each segment, when its length is no more than \eqn{vT}, the cost value
#' will be computed by performing an exact minimization of the loss function
#' over the parameter. When its length is greater than \eqn{vT}, the cost value
#' is approximated through SeGD. Therefore, this parameter induces an algorithm
#' that can be interpreted as an interpolation between dynamic programming with
#' SeGD (\eqn{v = 0}) and the vanilla PELT (\eqn{v = 1}).
#' The readers are referred to the references for more details.
#' @param warm_start If \code{TRUE}, the algorithm will use the estimated
#' parameters from the previous segment as the initial value for the
#' current segment. This parameter is only used for the \code{"glm"} families.
#' @param ... Other parameters for specific models.
#' \itemize{
#' \item \code{include.mean} is used to determine if a mean/intercept term
#' should be included in the ARIMA(\eqn{p}, \eqn{d}, \eqn{q}) or
#' GARCH(\eqn{p}, \eqn{q}) models.
#' \item \code{r.clock} is used to create an \code{RcppClock} object to record
#' the time spent in the C++ code. Default is an empty string. If set to
#' any non-empty string, an object with specified name will be created.
#' Usage: \code{library(RcppClock); plot(VARIABLE_NAME)}.
#' \item \code{r.progress} is used to control the progress bar. By default the
#' progress bar will be shown. To disable it, set \code{r.progress = FALSE}.
#' \item \code{p.response} is used to specify the number of response variables.
#' This parameter is especially useful for linear models with multivariate
#' responses.
#' }
#' @return A [fastcpd-class] object.
#' @description [fastcpd()] takes in formulas, data, families and extra
#' parameters and returns a [fastcpd-class] object.
#' @section Gallery:
#' <https://github.com/doccstat/fastcpd/tree/main/tests/testthat/examples>
#' @section References:
#' Xingchi Li, Xianyang Zhang (2024). ``fastcpd: Fast Change Point Detection
#' in R.'' _arXiv:2404.05933_, <https://arxiv.org/abs/2404.05933>.
#'
#' Xianyang Zhang, Trisha Dawn (2023). ``Sequential Gradient Descent and
#' Quasi-Newton's Method for Change-Point Analysis.'' In Ruiz, Francisco,
#' Dy, Jennifer, van de Meent, Jan-Willem (eds.), _Proceedings of The 26th
#' International Conference on Artificial Intelligence and Statistics_,
#' volume 206 series Proceedings of Machine Learning Research, 1129-1143.
#' @example tests/testthat/examples/fastcpd_1.R
#' @example tests/testthat/examples/fastcpd_2.R
#' @example tests/testthat/examples/fastcpd_3.txt
#' @example tests/testthat/examples/fastcpd_4.txt
#' @seealso [fastcpd.family] for the family-specific function;
#' [plot.fastcpd()] for plotting the results,
#' [summary.fastcpd()] for summarizing the results.
#'
#' @md
#' @importFrom fastglm fastglm
#' @importFrom glmnet glmnet cv.glmnet predict.glmnet
#' @importFrom methods show
#' @importFrom Rcpp evalCpp
#' @export
#' @useDynLib fastcpd, .registration = TRUE
fastcpd <- function(  # nolint: cyclomatic complexity
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
  trim = 0.02,
  momentum_coef = 0,
  multiple_epochs = function(x) 0,
  epsilon = 1e-10,
  order = c(0, 0, 0),
  p = ncol(data) - 1,
  cp_only = FALSE,
  vanilla_percentage = 0,
  warm_start = FALSE,
  ...
) {
  # Check the validity of the `family` parameter.
  check_family(
    family <- ifelse(is.null(family), "custom", tolower(family)),
    c(
      "lm",  # -> "gaussian"
      "binomial",  # -> "binomial"
      "poisson",  # -> "poisson"
      "lasso",  # -> "lasso"
      "mean",  # -> "mean"
      "variance",  # -> "variance"
      "meanvariance",  # -> "meanvariance"
      "arma",  # -> "arma"
      "ar",  # -> "gaussian"
      "var",  # -> "mgaussian"
      "arima",  # -> "custom"
      "garch",  # -> "custom"
      "custom"  # -> "custom"
    )
  )

  # Check the validity of the `cost` parameter.
  check_cost(cost, cost_gradient, cost_hessian, family)

  # Check the validity of the `cost_adjustment` parameter.
  if (is.null(cost_adjustment)) {
    cost_adjustment <- "BIC"
  }
  stopifnot(cost_adjustment %in% c("BIC", "MBIC", "MDL"))

  # The following code is adapted from the `lm` function from base R.
  match_formula <- match.call(expand.dots = FALSE)
  matched_formula <- match(c("formula", "data"), names(match_formula), 0L)
  match_formula <- match_formula[c(1L, matched_formula)]
  match_formula$drop.unused.levels <- TRUE
  match_formula[[1L]] <- quote(stats::model.frame)
  match_formula <- eval(match_formula, parent.frame())
  y <- stats::model.response(match_formula, "numeric")
  data_ <- cbind(y, stats::model.matrix(formula, data = data))

  if (family == "ar") {
    stopifnot("Data should be a univariate time series." = ncol(data_) == 1)
    stopifnot(check_ar_order(order))
  } else if (family == "var") {
    stopifnot(check_var_order(order))
  } else if (family == "garch") {
    stopifnot("Data should be a univariate time series." = ncol(data_) == 1)
    stopifnot(check_garch_order(order))
  } else if (family == "arima") {
    stopifnot("Data should be a univariate time series." = ncol(data_) == 1)
    stopifnot(check_arima_order(order))
  }

  # Check the parameters passed in the ellipsis.
  include_mean <- TRUE
  p_response <- get_p_response(family, y, data_)
  r_progress <- TRUE
  r_clock <- ""
  if (methods::hasArg("include.mean")) {
    include_mean <- eval.parent(match.call()[["include.mean"]])
  }
  if (methods::hasArg("p.response")) {
    p_response <- eval.parent(match.call()[["p.response"]])
  }
  if (methods::hasArg("r.progress")) {
    r_progress <- eval.parent(match.call()[["r.progress"]])
  }
  if (methods::hasArg("r.clock")) {
    r_clock <- eval.parent(match.call()[["r.clock"]])
  }

  # Define the cost functions for ARIMA and GARCH models.
  if (family == "arima") {
    cost <- function(data) {
      tryCatch(
        expr = -forecast::Arima(
          c(data), order = order, method = "ML", include.mean = include_mean
        )$loglik,
        error = function(e) 0
      )
    }
  } else if (family == "garch") {
    cost <- function(data) {
      tryCatch(
        expr = tseries::garch(data, order, trace = FALSE)$n.likeli,
        error = function(e) 0
      )
    }
  }

  # Get dimension of the data.
  d <- get_d(data_, family)

  # Get the number of paramters for the model to calculate the penalty.
  p <- get_p(data_, family, p_response, order, include_mean)

  # Assign families as "gaussian" for "lm" and "ar" or "mgaussian" for
  # "mlm" and "var".
  fastcpd_family <- get_fastcpd_family(family, p_response)

  # Estimate the variance / covariance matrix and pre-process the data for
  # mean, variance, meanvariance, ar and var models.
  sigma_data <- get_sigma_data(data_, family, order, p, p_response)
  sigma_ <- sigma_data$sigma
  data_ <- sigma_data$data

  # Use vanilla PELT for
  # "mean", "variance", "meanvariance", "arima", "garch", "mgaussian".
  vanilla_percentage <-
    get_vanilla_percentage(vanilla_percentage, cost, fastcpd_family)

  # Adjust the penalty for "lm".
  beta <- get_beta(beta, p, nrow(data_), fastcpd_family, sigma_)

  # No pruning for "lasso" and "mgaussian". Adjust the pruning coefficient for
  # "MBIC" and "MDL".
  pruning_coef <- get_pruning_coef(
    methods::hasArg("pruning_coef"),
    pruning_coef,
    cost_adjustment,
    fastcpd_family,
    p
  )

  result <- fastcpd_impl(
    data_, beta, cost_adjustment,  d, segment_count, trim, momentum_coef,
    multiple_epochs, fastcpd_family, epsilon, p, order, cost, cost_gradient,
    cost_hessian, cp_only, vanilla_percentage, warm_start, lower, upper,
    line_search, sigma_, p_response, pruning_coef, r_clock, r_progress
  )

  raw_cp_set <- c(result$raw_cp_set)
  cp_set <- c(result$cp_set)

  if (family %in% c("ar", "var") && length(order) == 1) {
    raw_cp_set <- raw_cp_set + order
    cp_set <- cp_set + order
  }

  if (vanilla_percentage == 1) {
    thetas <- data.frame(matrix(NA, 0, 0))
  } else {
    thetas <- data.frame(result$thetas)
    if (ncol(thetas) > 0) {
      names(thetas) <- paste0("segment ", seq_len(ncol(thetas)))
    }
  }

  if (is.null(result$cost_values)) {
    result$cost_values <- numeric(0)
  }

  if (is.null(result$residual)) {
    result$residual <- numeric(0)
  }

  residuals <- matrix(result$residual)

  if (!cp_only) {
    tryCatch(
      expr = if (family == "ar") {
        residuals <- matrix(c(rep(NA, p), residuals))
      } else if (family == "arima") {
        residuals <- matrix(NA, nrow(data))
        segments <- c(0, cp_set, nrow(data))
        for (segments_i in seq_len(length(segments) - 1)) {
          segments_start <- segments[segments_i] + 1
          segments_end <- segments[segments_i + 1]
          residuals[segments_start:segments_end] <- forecast::Arima(
            c(data[segments_start:segments_end, 1]),
            order = order,
            method = "ML",
            include.mean = include_mean
          )$residuals
        }
      } else if (family == "mean") {
        residuals <- matrix(NA, nrow(data), ncol(data))
        segments <- c(0, cp_set, nrow(data))
        for (segments_i in seq_len(length(segments) - 1)) {
          segments_start <- segments[segments_i] + 1
          segments_end <- segments[segments_i + 1]
          sgmt_index <- segments_start:segments_end
          residuals[sgmt_index, ] <- as.matrix(
            data[sgmt_index, ] - colMeans(data[sgmt_index, , drop = FALSE])
          )
        }
      } else if (family == "garch") {
        residuals <- matrix(NA, nrow(data))
        segments <- c(0, cp_set, nrow(data))
        for (segments_i in seq_len(length(segments) - 1)) {
          segments_start <- segments[segments_i] + 1
          segments_end <- segments[segments_i + 1]
          residuals[segments_start:segments_end] <- tseries::garch(
            data[segments_start:segments_end, 1],
            order,
            trace = FALSE
          )$residuals
        }
      },
      error = function(e) message("Residual calculation failed.")
    )
  }

  methods::new(
    Class = "fastcpd",
    call = match.call(),
    data = data.frame(data),
    order = order,
    family = family,
    cp_set = cp_set,
    cost_values = c(result$cost_values),
    residuals = residuals,
    thetas = thetas,
    cp_only = cp_only
  )
}
