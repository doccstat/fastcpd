#' @title Find change points efficiently
#'
#' @description \code{fastcpd} takes in formulas, data, families and extra
#' parameters and returns a \code{fastcpd} object.
#'
#' @example tests/testthat/examples/fastcpd.txt
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
#' @param family Family of the model. Can be \code{"lm"}, \code{"binomial"},
#'   \code{"poisson"}, \code{"lasso"}, \code{"custom"}, \code{"ar"},
#'   \code{"var"}, \code{"ma"}, \code{"arima"}, \code{"garch"} or
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
#'   \code{p = ncol(data) - 1}. This parameter is superseded by `order` in the
#'   case of time series models: "ar", "var", "arima".
#' @param order Order of the AR(p), VAR(p) or ARIMA(p, d, q) model.
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
#' @param ... Parameters specifically used for time series models. As of the
#'   current implementation, only `include.mean` will not be ignored and used
#'   in the ARIMA or GARCH model.
#'
#' @return A class \code{fastcpd} object.
#'
#' @export
#' @importFrom DescTools Winsorize
#' @importFrom fastglm fastglm
#' @importFrom glmnet glmnet cv.glmnet predict.glmnet
#' @importFrom methods show
#' @importFrom Rcpp evalCpp
#' @useDynLib fastcpd, .registration = TRUE
fastcpd <- function(  # nolint: cyclomatic complexity
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
  p = ncol(data) - 1,
  order = c(0, 0, 0),
  cost = NULL,
  cost_gradient = NULL,
  cost_hessian = NULL,
  cp_only = FALSE,
  vanilla_percentage = 0,
  warm_start = FALSE,
  ...
) {
  family <- ifelse(is.null(family), "custom", tolower(family))

  # Vanilla is not a `fastcpd` family and can not be set manually by the user.
  # `vanilla` is used to distinguish the cost function parameters in the
  # implementation.
  stopifnot(
    check_family(
      family,
      c(
        "lm", "binomial", "poisson", "lasso",
        "ar", "var", "ma", "arima", "garch", "custom"
      )
    )
  )

  stopifnot(check_cost(cost, cost_gradient, cost_hessian, family))

  # The following code is adapted from the `lm` function from base R.
  match_formula <- match.call(expand.dots = FALSE)
  matched_formula <- match(c("formula", "data"), names(match_formula), 0L)
  match_formula <- match_formula[c(1L, matched_formula)]
  match_formula$drop.unused.levels <- TRUE
  match_formula[[1L]] <- quote(stats::model.frame)
  match_formula <- eval(match_formula, parent.frame())
  y <- stats::model.response(match_formula, "any")
  data <- cbind(y, stats::model.matrix(formula, data = data))

  fastcpd_family <- NULL
  fastcpd_data <- NULL

  # If a cost function provided has an explicit solution, i.e. does not depend
  # on the parameters, e.g., mean change, then the percentage of vanilla PELT
  # is set to be 1.
  if (!is.null(cost) && length(formals(cost)) == 1) {
    family <- "custom"
    vanilla_percentage <- 1
  }

  # Pre-process the data for the time series models.
  if (family == "ar") {
    # Check the validity of the parameters for AR(p) model.
    stopifnot("Data should be a univariate time series." = ncol(data) == 1)
    stopifnot(check_ar_order(order))

    if (length(order) == 3) {
      family <- "arima"
    } else {
      # Preprocess the data for AR(p) model to be used in linear regression.
      p <- order
      fastcpd_family <- "gaussian"
      y <- data[p + seq_len(nrow(data) - p), ]
      x <- matrix(NA, nrow(data) - p, p)
      for (p_i in seq_len(p)) {
        x[, p_i] <- data[(p - p_i) + seq_len(nrow(data) - p), ]
      }
      fastcpd_data <- cbind(y, x)
    }
  } else if (family == "var") {
    stopifnot(check_var_order(order))
    fastcpd_family <- "custom"

    # Preprocess the data for VAR(p) model to be used in linear regression.
    vanilla_percentage <- 1
    y <- data[order + seq_len(nrow(data) - order), ]
    x <- matrix(NA, nrow(data) - order, order * ncol(data))
    for (p_i in seq_len(order)) {
      x[, (p_i - 1) * ncol(data) + seq_len(ncol(data))] <-
        data[(order - p_i) + seq_len(nrow(data) - order), ]
    }
    fastcpd_data <- cbind(y, x)
    lm_x_col <- order * ncol(data)
    cost <- function(data) {
      x <- data[, (ncol(data) - lm_x_col + 1):ncol(data)]
      y <- data[, 1:(ncol(data) - lm_x_col)]

      if (nrow(data) <= lm_x_col + 1) {
        x_t_x <- diag(lm_x_col)
      } else {
        x_t_x <- crossprod(x)
      }

      # TODO(doccstat): Verify the correctness of the cost function for VAR(p).
      norm(y - x %*% solve(x_t_x, t(x)) %*% y, type = "F")^2 / 2
    }
    p <- order^2 * ncol(data)
  } else if (family == "ma") {
    stopifnot("Data should be a univariate time series." = ncol(data) == 1)
    stopifnot(check_ma_order(order))
    fastcpd_family <- "custom"
    vanilla_percentage <- 1

    if (length(order) == 1) {
      order <- c(0, 0, order)
    }

    # By default in time series models, `include.mean` is set to be `TRUE`.
    include_mean <- TRUE
    if (methods::hasArg("include.mean")) {
      include_mean <- eval.parent(match.call()[["include.mean"]])
    }
    cost <- function(data) {
      tryCatch(
        expr = -forecast::Arima(
          c(data), order = order, method = "ML", include.mean = include_mean
        )$loglik,
        error = function(e) 0
      )
    }
    p <- sum(order)
  }
  if (family == "arima") {
    stopifnot("Data should be a univariate time series." = ncol(data) == 1)
    stopifnot(check_arima_order(order))
    fastcpd_family <- "custom"
    vanilla_percentage <- 1

    # By default in time series models, `include.mean` is set to be `TRUE`.
    include_mean <- TRUE
    if (methods::hasArg("include.mean")) {
      include_mean <- eval.parent(match.call()[["include.mean"]])
    }
    cost <- function(data) {
      tryCatch(
        expr = -forecast::Arima(
          c(data), order = order, method = "ML", include.mean = include_mean
        )$loglik,
        error = function(e) 0
      )
    }
    p <- sum(order[-2])
  } else if (family == "garch") {
    stopifnot("Data should be a univariate time series." = ncol(data) == 1)
    stopifnot(check_garch_order(order))
    fastcpd_family <- "custom"
    vanilla_percentage <- 1

    # By default in time series models, `include.mean` is set to be `TRUE`.
    include_mean <- TRUE
    if (methods::hasArg("include.mean")) {
      include_mean <- eval.parent(match.call()[["include.mean"]])
    }
    garch_formula <-
      stats::as.formula(paste0("~ garch(", order[1], ",", order[2], ")"))
    cost <- function(data) {
      tryCatch(
        expr = fGarch::garchFit(
          formula = garch_formula,
          data = data,
          include.mean = include_mean,
          trace = FALSE
        )@fit$value,
        error = function(e) 0
      )
    }
    p <- 1 + sum(order)
  }

  if (family == "lm") {
    fastcpd_family <- "gaussian"
  }

  if (is.null(fastcpd_family)) {
    fastcpd_family <- family
  }

  if (is.null(fastcpd_data)) {
    fastcpd_data <- data
  }

  if (is.null(beta)) {
    # Use the `beta` value obtained from BIC.
    beta <- (p + 1) * log(nrow(fastcpd_data)) / 2

    # TODO(doccstat): Variance estimation for VAR(p).
    # For linear regression models, an estimate of the variance is needed in the
    # cost function. The variance estimation is only for "lm" family with no
    # `beta` provided. Only estimate the variance for Gaussian family when
    # `beta` is null.
    if (family == "lm" || fastcpd_family == "gaussian") {
      # Estimate the variance for each block and then take the average.
      n <- nrow(fastcpd_data)
      block_size <- 5
      variance_estimation <- rep(NA, n - block_size)
      for (i in 1:(n - block_size)) {
        block_index <- seq_len(block_size) + i - 1
        block_index_lagged <- seq_len(block_size) + i

        y_block <- fastcpd_data[block_index, 1]
        x_block <- fastcpd_data[block_index, -1, drop = FALSE]

        y_block_lagged <- fastcpd_data[block_index_lagged, 1]
        x_block_lagged <- fastcpd_data[block_index_lagged, -1, drop = FALSE]

        x_t_x <- crossprod(x_block)
        x_t_x_lagged <- crossprod(x_block_lagged)

        block_slope <- solve(crossprod(x_block), crossprod(x_block, y_block))
        block_lagged_slope <- solve(
          crossprod(x_block_lagged), crossprod(x_block_lagged, y_block_lagged)
        )
        x_t_x_inv <- solve(x_t_x)
        x_t_x_inv_lagged <- solve(x_t_x_lagged)
        inv_product <- x_t_x_inv %*% x_t_x_inv_lagged
        cross_term_x <- crossprod(
          x_block[-1, , drop = FALSE],
          x_block_lagged[-block_size, , drop = FALSE]
        )
        cross_term <- inv_product %*% cross_term_x
        delta_numerator <- crossprod(block_slope - block_lagged_slope)
        delta_denominator <-
          sum(diag(x_t_x_inv + x_t_x_inv_lagged - 2 * cross_term))
        variance_estimation[i] <- delta_numerator / delta_denominator
      }

      beta <- beta * mean(variance_estimation)
    }
  }

  result <- fastcpd_impl(
    fastcpd_data, beta, segment_count, trim, momentum_coef, k, fastcpd_family,
    epsilon, min_prob, winsorise_minval, winsorise_maxval, p,
    cost, cost_gradient, cost_hessian, cp_only, vanilla_percentage, warm_start
  )

  cp_set <- c(result$cp_set)

  if (family %in% c("ar", "var") && length(order) == 1) {
    cp_set <- cp_set + order
  }

  if (vanilla_percentage == 1) {
    thetas <- data.frame(matrix(NA, 0, 0))
  } else {
    thetas <- data.frame(result$thetas)
    names(thetas) <- paste0("segment ", seq_len(ncol(thetas)))
  }

  if (is.null(result$cost_values)) {
    result$cost_values <- numeric(0)
  }

  if (is.null(result$residual)) {
    result$residual <- numeric(0)
  }

  residuals <- c(result$residual)

  if (family == "ar") {
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
      residuals[segments_start:segments_end] <- fGarch::garchFit(
        formula = garch_formula,
        data = data[segments_start:segments_end, 1],
        include.mean = include_mean,
        trace = FALSE
      )@residuals
    }
  }

  methods::new(
    Class = "fastcpd",
    call = match.call(),
    data = data.frame(data),
    family = family,
    cp_set = cp_set,
    cost_values = c(result$cost_values),
    residuals = residuals,
    thetas = thetas,
    cp_only = cp_only
  )
}
