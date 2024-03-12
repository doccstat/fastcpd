#' @title Find change points efficiently
#'
#' @description \code{fastcpd} takes in formulas, data, families and extra
#' parameters and returns a \code{fastcpd} object.
#'
#' @example tests/testthat/examples/fastcpd_1.R
#' @example tests/testthat/examples/fastcpd_2.R
#' @example tests/testthat/examples/fastcpd_3.txt
#' @example tests/testthat/examples/fastcpd_4.R
#'
#' @md
#' @section Gallery:
#'   <https://fastcpd.xingchi.li/articles/gallery.html>
#' @section References:
#'   Zhang X, Dawn T (2023). ``Sequential Gradient Descent and Quasi-Newton's
#'   Method for Change-Point Analysis.'' In Ruiz, Francisco, Dy, Jennifer,
#'   van de Meent, Jan-Willem (eds.), _Proceedings of The 26th International
#'   Conference on Artificial Intelligence and Statistics_, volume 206 series
#'   Proceedings of Machine Learning Research, 1129-1143.
#'   <https://proceedings.mlr.press/v206/zhang23b.html>.
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
#' @param beta Penalty criterion for the number of change points.
#'   Take a string value of \code{"MBIC"}, \code{"BIC"}, \code{"MDL"}
#'   or a numeric value. If a numeric value is provided, the value will be
#'   used as the penalty criterion for the number of change points in the
#'   traditional Bayesian Information Criterion. By default, Modified BIC
#'   criterion is used to obtain a proper value, i.e.,
#'   \code{beta = (p + 2) * log(nrow(data)) / 2} with adjusted cost
#'   function.
#' @param cost_adjustment Cost adjustment criterion for the number of change
#'   points. Can be \code{"BIC"}, \code{"MBIC"}, \code{"MDL"} or NULL.
#'   By default, the cost adjustment criterion is set to be \code{"MBIC"}.
#'   MBIC and MDL modifies the cost function by adding a small negative
#'   term to the cost function. MDL then transforms the cost function to
#'   log2 based. BIC or NULL does not modify the cost function.
#' @param family Family of the model. Can be \code{"lm"}, \code{"binomial"},
#'   \code{"poisson"}, \code{"lasso"}, \code{"custom"}, \code{"ar"},
#'   \code{"var"}, \code{"ma"}, \code{"arima"}, \code{"garch"} or
#'   \code{NULL}. For simplicity, user can also omit this parameter,
#'   indicating that they will be using their own cost functions. Omitting the
#'   parameter is the same as specifying the parameter to be \code{"custom"}
#'   or \code{NULL}, in which case, users must specify the cost function, with
#'   optional gradient and corresponding Hessian matrix functions.
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
#' @param line_search If a vector of numeric values are provided, line
#'   search will be performed to find the optimal step size for each update.
#' @param lower Lower bound for the parameters. Used to specify the
#'   domain of the parameter after each gradient descent step. If not specified,
#'   the lower bound will be set to be \code{-Inf} for all parameters.
#' @param upper Upper bound for the parameters. Used to specify the
#'   domain of the parameter after each gradient descent step. If not specified,
#'   the upper bound will be set to be \code{Inf} for all parameters.
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
#' @param multiple_epochs Function on number of epochs in SGD.
#'   \code{multiple_epochs} should be a function taking only a parameter
#'   \code{segment_length} meaning the current number of data
#'   points considered since last segmentaion. The return value of the
#'   function should be an integer indicating how many epochs should be
#'   performed apart from the default update. By default the function returns
#'   0, meaning no multiple epochs will be used to update the parameters.
#'   Example usage:
#'   ```r
#'     multiple_epochs = function(segment_length) {
#'       if (segment_length < n / segment_count / 4 * 1) 3
#'       else if (segment_length < n / segment_count / 4 * 2) 2
#'       else if (segment_length < n / segment_count / 4 * 3) 1
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
#' @param epsilon Epsilon to avoid numerical issues. Only used for the Hessian
#'   computation in Logistic Regression and Poisson Regression.
#' @param order Order of the AR(p), VAR(p) or ARIMA(p, d, q) model.
#' @param p Number of covariates in the model. If not specified, the number of
#'   covariates will be inferred from the data, i.e.,
#'   \code{p = ncol(data) - 1}. This parameter is superseded by `order` in the
#'   case of time series models: "ar", "var", "arima".
#' @param pruning If \code{TRUE}, the algorithm will perform pruning on the
#'   change points. By default, the value is set to be \code{TRUE}. Pruning
#'   should be set to be \code{FALSE} if the pruning condition is not
#'   satisfied.
#' @param cp_only If \code{TRUE}, only the change points are returned.
#'   Otherwise, the cost function values together with the estimated
#'   parameters for each segment are also returned. By default the value is
#'   set to be \code{FALSE} so that `plot` can be used to visualize the
#'   results for a built-in model. \code{cp_only} has some performance impact
#'   on the algorithm, since the cost values and estimated parameters for each
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
#' @param warm_start If \code{TRUE}, the algorithm will use the estimated
#'   parameters from the previous segment as the initial value for the
#'   current segment. This parameter is only used for \code{"glm"} families.
#' @param ... Parameters specifically used for time series models. The following
#'   parameters will not be ignored:
#'
#'   - \code{include.mean} is used for ARIMA and GARCH modes.
#'   - \code{r.progress} is used to control the progress bar. By default the
#'     progress bar will be shown. To disable it, set \code{r.progress = FALSE}.
#'   - \code{p.response} is used to specify the number of response variables.
#'
#' @return A class \code{fastcpd} object.
#'
#' @export
#' @importFrom fastglm fastglm
#' @importFrom glmnet glmnet cv.glmnet predict.glmnet
#' @importFrom methods show
#' @importFrom Rcpp evalCpp
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
  segment_count = 10,
  trim = 0.02,
  momentum_coef = 0,
  multiple_epochs = function(x) 0,
  epsilon = 1e-10,
  order = c(0, 0, 0),
  p = ncol(data) - 1,
  pruning = TRUE,
  cp_only = FALSE,
  vanilla_percentage = 0,
  warm_start = FALSE,
  ...
) {
  # Check the validity of the `family` parameter.
  family <- ifelse(is.null(family), "custom", tolower(family))
  stopifnot(
    check_family(
      family,
      c(
        "lm", "binomial", "poisson", "lasso", "mlasso", "mean", "variance",
        "meanvariance", "mv", "arma", "ar", "var", "ma", "arima", "garch",
        "custom"
      )
    )
  )

  # Check the validity of the `cost` parameter.
  stopifnot(check_cost(cost, cost_gradient, cost_hessian, family))

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

  # Check the parameters passed in the ellipsis.
  include_mean <- TRUE
  p_response <- get_p_response(family, y, data)
  r_progress <- TRUE
  if (methods::hasArg("include.mean")) {
    include_mean <- eval.parent(match.call()[["include.mean"]])
  }
  if (methods::hasArg("p.response")) {
    p_response <- eval.parent(match.call()[["p.response"]])
  }
  if (methods::hasArg("r.progress")) {
    r_progress <- eval.parent(match.call()[["r.progress"]])
  }

  p <- get_p(data_, family, p_response, order, include_mean)

  if (family == "ar") {
    stopifnot("Data should be a univariate time series." = ncol(data_) == 1)
    stopifnot(check_ar_order(order))
  } else if (family == "var") {
    stopifnot(check_var_order(order))
  } else if (family == "ma") {
    stopifnot("Data should be a univariate time series." = ncol(data_) == 1)
    stopifnot(check_ma_order(order))
  } else if (family == "garch") {
    stopifnot("Data should be a univariate time series." = ncol(data_) == 1)
    stopifnot(check_garch_order(order))
  } else if (family == "arima") {
    stopifnot("Data should be a univariate time series." = ncol(data_) == 1)
    stopifnot(check_arima_order(order))
  }

  if (family == "ar") {
    y <- data_[p + seq_len(nrow(data_) - p), ]
    x <- matrix(NA, nrow(data_) - p, p)
    for (p_i in seq_len(p)) {
      x[, p_i] <- data_[(p - p_i) + seq_len(nrow(data_) - p), ]
    }
    data_ <- cbind(y, x)
  } else if (family == "var") {
    y <- data_[order + seq_len(nrow(data_) - order), ]
    x <- matrix(NA, nrow(data_) - order, order * ncol(data_))
    for (p_i in seq_len(order)) {
      x[, (p_i - 1) * ncol(data_) + seq_len(ncol(data_))] <-
        data_[(order - p_i) + seq_len(nrow(data_) - order), ]
    }
    data_ <- cbind(y, x)
  } else if (family == "ma") {
    # TODO(doccstat): Deprecate MA model.
    family <- "arima"
    order <- c(rep(0, 3 - length(order)), order)
  } else if (family == "garch") {
    cost <- function(data) {
      tryCatch(
        expr = tseries::garch(data, order, trace = FALSE)$n.likeli,
        error = function(e) 0
      )
    }
  }

  if (family == "arima") {
    cost <- function(data) {
      tryCatch(
        expr = -forecast::Arima(
          c(data), order = order, method = "ML", include.mean = include_mean
        )$loglik,
        error = function(e) 0
      )
    }
  }

  fastcpd_family <- get_fastcpd_family(family, p_response)
  sigma_ <- get_variance_estimation(data_, family, p_response)
  vanilla_percentage <-
    get_vanilla_percentage(vanilla_percentage, cost, fastcpd_family)
  beta <- get_beta(beta, p, nrow(data_), fastcpd_family, sigma_)

  result <- fastcpd_impl(
    data_, beta, cost_adjustment, segment_count, trim, momentum_coef,
    multiple_epochs, fastcpd_family, epsilon, p, pruning, order, cost,
    cost_gradient, cost_hessian, cp_only, vanilla_percentage, warm_start, lower,
    upper, line_search, sigma_, p_response, r_progress
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

  raw_residuals <- c(result$residual)

  if (family == "mean" && p == 1) {
    segments <- c(0, raw_cp_set, nrow(data))
    for (segments_i in seq_len(length(segments) - 1)) {
      segments_start <- segments[segments_i] + 1
      segments_end <- segments[segments_i + 1]
      segment_index <- segments_start:segments_end
      raw_residuals[segment_index] <-
        data[segment_index, 1] - mean(data[segment_index, 1])
    }
  }

  if (!cp_only) {
    tryCatch(
      expr = if (family == "ar") {
        raw_residuals <- c(rep(NA, p), raw_residuals)
      } else if (family == "ma" || family == "arima") {
        raw_residuals <- rep(NA, nrow(data))
        segments <- c(0, raw_cp_set, nrow(data))
        for (segments_i in seq_len(length(segments) - 1)) {
          segments_start <- segments[segments_i] + 1
          segments_end <- segments[segments_i + 1]
          raw_residuals[segments_start:segments_end] <- forecast::Arima(
            c(data[segments_start:segments_end, 1]),
            order = order,
            method = "ML",
            include.mean = include_mean
          )$residuals
        }
      } else if (family == "garch") {
        raw_residuals <- rep(NA, nrow(data))
        segments <- c(0, raw_cp_set, nrow(data))
        for (segments_i in seq_len(length(segments) - 1)) {
          segments_start <- segments[segments_i] + 1
          segments_end <- segments[segments_i + 1]
          raw_residuals[segments_start:segments_end] <- tseries::garch(
            data[segments_start:segments_end, 1],
            order,
            trace = FALSE
          )$residuals
        }
      },
      error = function(e) message("Residual calculation failed.")
    )
  }

  residuals <- c(result$residual)

  if (family == "mean" && p == 1) {
    segments <- c(0, cp_set, nrow(data))
    for (segments_i in seq_len(length(segments) - 1)) {
      segments_start <- segments[segments_i] + 1
      segments_end <- segments[segments_i + 1]
      segment_index <- segments_start:segments_end
      residuals[segment_index] <-
        data[segment_index, 1] - mean(data[segment_index, 1])
    }
  }

  if (!cp_only) {
    tryCatch(
      expr = if (family == "ar") {
        residuals <- c(rep(NA, p), residuals)
      } else if (family == "ma" || family == "arima") {
        residuals <- rep(NA, nrow(data))
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
      } else if (family == "garch") {
        residuals <- rep(NA, nrow(data))
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
