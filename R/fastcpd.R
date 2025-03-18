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
#' @param variance_estimation An estimate of the variance / covariance matrix
#' for the data. If not specified, the variance / covariance matrix will be
#' estimated using the data.
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
  variance_estimation = NULL,
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
      "garch",  # -> "garch"
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
  if (methods::hasArg("include.mean")) {
    include_mean <- eval.parent(match.call()[["include.mean"]])
  }
  if (methods::hasArg("p.response")) {
    p_response <- eval.parent(match.call()[["p.response"]])
  }
  if (methods::hasArg("r.progress")) {
    r_progress <- eval.parent(match.call()[["r.progress"]])
  }

  if (family %in% c("binomial", "poisson", "lasso")) {
    fastcpd_family <- family
  } else if (family == "mean") {
    fastcpd_family <- family
    vanilla_percentage <- 1
    p <- ncol(data_)
  } else if (family == "variance") {
    fastcpd_family <- family
    vanilla_percentage <- 1
    p <- ncol(data_)^2
  } else if (family == "meanvariance") {
    fastcpd_family <- family
    vanilla_percentage <- 1
    p <- ncol(data_)^2 + ncol(data_)
  } else if (family == "garch") {
    p <- sum(order) + 1
    fastcpd_family <- family
    vanilla_percentage <- 1
  } else if (family == "lm" && p_response == 1) {
    fastcpd_family <- "gaussian"
  } else if (family == "ar") {
    p <- order
    fastcpd_family <- "gaussian"
    y <- data_[p + seq_len(nrow(data_) - p), ]
    x <- matrix(NA, nrow(data_) - p, p)
    for (p_i in seq_len(p)) {
      x[, p_i] <- data_[(p - p_i) + seq_len(nrow(data_) - p), ]
    }
    data_ <- cbind(y, x)
  } else if (family == "lm" && p_response > 1) {
    p <- (ncol(data_) - p_response) * p_response
    fastcpd_family <- "mgaussian"
    vanilla_percentage <- 1
  } else if (family == "var") {
    p <- order * p_response^2
    fastcpd_family <- "mgaussian"
    vanilla_percentage <- 1
    y <- data_[order + seq_len(nrow(data_) - order), ]
    x <- matrix(NA, nrow(data_) - order, order * ncol(data_))
    for (p_i in seq_len(order)) {
      x[, (p_i - 1) * ncol(data_) + seq_len(ncol(data_))] <-
        data_[(order - p_i) + seq_len(nrow(data_) - order), ]
    }
    data_ <- cbind(y, x)
  } else if (family == "arma" && order[1] == 0) {
    p <- sum(order) + 1
    fastcpd_family <- "ma"
  } else if (family == "arma" && order[1] != 0) {
    p <- sum(order) + 1
    fastcpd_family <- family
  } else if (family == "arima") {
    cost <- function(data) {
      tryCatch(
        expr = -stats::arima(
          c(data), order = order, method = "ML", include.mean = include_mean
        )$loglik,
        error = function(e) 0
      )
    }
    p <- sum(order[-2]) + 1 + include_mean
    fastcpd_family <- "custom"
    if (!is.null(cost) && length(formals(cost)) == 1) {
      vanilla_percentage <- 1
    }
  } else {
    if (!methods::hasArg("p")) {
      p <- ncol(data_) - 1
    }
    fastcpd_family <- "custom"
    if (!is.null(cost) && length(formals(cost)) == 1) {
      vanilla_percentage <- 1
    }
  }

  cost_pelt <- NULL
  cost_sen <- NULL
  if (length(formals(cost)) == 1) {
    cost_pelt <- cost
  } else {
    cost_sen <- cost
  }

  sigma_ <- if (!is.null(variance_estimation)) {
    as.matrix(variance_estimation)
  } else if (family == "mean") {
    variance.mean(data_)
  } else if (family == "var" || family == "lm" && p_response > 1) {
    as.matrix(Matrix::nearPD(variance.lm(data_, p_response))$mat)
  } else if (family == "lm" || family == "ar") {
    as.matrix(variance.lm(data_))
  } else {
    diag(1)
  }

  if (rcond(sigma_) < 1e-10) {
    sigma_ <- diag(1e-10, nrow(sigma_))
  }

  if (is.character(beta)) {
    if (!(beta %in% c("BIC", "MBIC", "MDL"))) {
      stop("Invalid beta selection criterion provided.")
    }

    beta <- switch(
      beta,
      BIC = (p + 1) * log(nrow(data_)) / 2,
      MBIC = (p + 2) * log(nrow(data_)) / 2,
      MDL = (p + 2) * log2(nrow(data_)) / 2
    )

    # For linear regression models, an estimate of the variance is needed in the
    # cost function. The variance estimation is only for "lm" family with no
    # `beta` provided. Only estimate the variance for Gaussian family when
    # `beta` is null.
    if (fastcpd_family == "gaussian") {
      beta <- beta * c(sigma_)
    }
  }

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
    data_, beta, cost_adjustment, segment_count, trim, momentum_coef,
    multiple_epochs, fastcpd_family, epsilon, p, order, cost_pelt, cost_sen,
    cost_gradient, cost_hessian, cp_only, vanilla_percentage, warm_start,
    lower, upper, line_search, sigma_, p_response, pruning_coef, r_progress
  )

  raw_cp_set <- c(result$raw_cp_set)
  cp_set <- c(result$cp_set)

  if (family %in% c("ar", "var") && length(order) == 1) {
    raw_cp_set <- raw_cp_set + order
    cp_set <- cp_set + order
  }

  thetas <- data.frame(result$thetas)
  if (ncol(thetas) > 0) {
    names(thetas) <- paste0("segment ", seq_len(ncol(thetas)))
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
          residuals[segments_start:segments_end] <- stats::arima(
            c(data[segments_start:segments_end, 1]),
            order = order,
            method = "ML",
            include.mean = include_mean
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

#' @name fastcpd_family
#' @aliases fastcpd.family
#' @title Wrapper functions for fastcpd
#' @description Wrapper functions for fastcpd to find change points in various
#' models.
#' @seealso [fastcpd.mean()], [fastcpd.variance()], [fastcpd.mv()],
#' [fastcpd.meanvariance()] for basic statistics change models;
#' [fastcpd.lm()], [fastcpd.binomial()], [fastcpd.poisson()],
#' [fastcpd.lasso()] for regression coefficients change models;
#' [fastcpd.ar()], [fastcpd.var()], [fastcpd.arima()], [fastcpd.arma()],
#' [fastcpd.garch()] for change in time series models.
#'
#' @md
#' @keywords internal
NULL

#' @title Find change points efficiently in AR(\eqn{p}) models
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param order A positive integer specifying the order of the AR model.
#' @param ... Other arguments passed to [fastcpd()], for example,
#' \code{segment_count}. One special argument can be passed here is
#' \code{include.mean}, which is a logical value indicating whether the
#' mean should be included in the model. The default value is \code{TRUE}.
#' @return A [fastcpd-class] object.
#' @description [fastcpd_ar()] and [fastcpd.ar()] are
#' wrapper functions of [fastcpd()] to find change points in
#' AR(\eqn{p}) models. The function is similar to [fastcpd()] except that
#' the data is by default a one-column matrix or univariate vector
#' and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_ar_1.R
#' @example tests/testthat/examples/fastcpd_ar_2.R
#' @seealso [fastcpd()]
#'
#' @md
#' @rdname fastcpd_ar
#' @export
fastcpd_ar <- function(data, order = 0, ...) {
  result <- fastcpd.ts(c(data), "ar", order, ...)
  result@call <- match.call()
  result
}

#' @rdname fastcpd_ar
#' @export
fastcpd.ar <- fastcpd_ar  # nolint: Conventional R function style

#' @title Find change points efficiently in
#' ARIMA(\eqn{p}, \eqn{d}, \eqn{q}) models
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param order A vector of length three specifying the order of the ARIMA
#' model.
#' @param ... Other arguments passed to [fastcpd()], for example,
#' \code{segment_count}. One special argument can be passed here is
#' \code{include.mean}, which is a logical value indicating whether the
#' mean should be included in the model. The default value is \code{TRUE}.
#' @return A [fastcpd-class] object.
#' @description [fastcpd_arima()] and [fastcpd.arima()] are
#' wrapper functions of [fastcpd()] to find change points in
#' ARIMA(\eqn{p}, \eqn{d}, \eqn{q}) models.
#' The function is similar to [fastcpd()]
#' except that the data is by default a one-column matrix or univariate vector
#' and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_arima.txt
#' @seealso [fastcpd()]
#'
#' @md
#' @rdname fastcpd_arima
#' @export
fastcpd_arima <- function(data, order = 0, ...) {
  result <- fastcpd.ts(c(data), "arima", order, ...)
  result@call <- match.call()
  result
}

#' @rdname fastcpd_arima
#' @export
fastcpd.arima <- fastcpd_arima  # nolint: Conventional R function style

#' @title Find change points efficiently in ARMA(\eqn{p}, \eqn{q}) models
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param order A vector of length two specifying the order of the ARMA
#' model.
#' @param ... Other arguments passed to [fastcpd()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [fastcpd_arma()] and [fastcpd.arma()] are
#' wrapper functions of [fastcpd()] to find change points in
#' ARMA(\eqn{p}, \eqn{q}) models. The function is similar to [fastcpd()]
#' except that the data is by default a one-column matrix or univariate vector
#' and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_arma.txt
#' @seealso [fastcpd()]
#'
#' @md
#' @rdname fastcpd_arma
#' @export
fastcpd_arma <- function(data, order = c(0, 0), ...) {
  result <- fastcpd.ts(c(data), "arma", order, ...)
  result@call <- match.call()
  result
}

#' @rdname fastcpd_arma
#' @export
fastcpd.arma <- fastcpd_arma  # nolint: Conventional R function style

#' @title Find change points efficiently in logistic regression models
#' @param data A matrix or a data frame with the response variable as the first
#' column.
#' @param ... Other arguments passed to [fastcpd()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [fastcpd_binomial()] and [fastcpd.binomial()] are
#' wrapper functions of [fastcpd()] to find change points in
#' logistic regression models. The function is similar to [fastcpd()]
#' except that the data is by default a matrix or data frame with the response
#' variable as the first column and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_binomial.R
#' @seealso [fastcpd()]
#'
#' @md
#' @rdname fastcpd_binomial
#' @export
fastcpd_binomial <- function(data, ...) {
  result <- fastcpd(
    data = data.frame(y = data[, 1], x = data[, -1]), family = "binomial", ...
  )
  result@call <- match.call()
  result
}

#' @rdname fastcpd_binomial
#' @export
fastcpd.binomial <- fastcpd_binomial  # nolint: Conventional R function style

#' @title Find change points efficiently in GARCH(\eqn{p}, \eqn{q}) models
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param order A positive integer vector of length two specifying the order of
#' the GARCH model.
#' @param ... Other arguments passed to [fastcpd()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [fastcpd_garch()] and [fastcpd.garch()] are
#' wrapper functions of [fastcpd()] to find change points in
#' GARCH(\eqn{p}, \eqn{q}) models. The function is similar to [fastcpd()]
#' except that the data is by default a one-column matrix or univariate vector
#' and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_garch.txt
#' @seealso [fastcpd()]
#'
#' @md
#' @rdname fastcpd_garch
#' @export
fastcpd_garch <- function(data, order = c(0, 0), ...) {
  result <- fastcpd.ts(c(data), "garch", order, ...)
  result@call <- match.call()
  result
}

#' @rdname fastcpd_garch
#' @export
fastcpd.garch <- fastcpd_garch  # nolint: Conventional R function style

#' @title Find change points efficiently in penalized linear regression models
#' @param data A matrix or a data frame with the response variable as the first
#' column.
#' @param ... Other arguments passed to [fastcpd()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [fastcpd_lasso()] and [fastcpd.lasso()] are wrapper
#' functions of [fastcpd()] to find change points in penalized
#' linear regression models. The function is similar to [fastcpd()]
#' except that the data is by default a matrix or data frame with the response
#' variable as the first column and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_lasso.txt
#' @seealso [fastcpd()]
#'
#' @md
#' @rdname fastcpd_lasso
#' @export
fastcpd_lasso <- function(data, ...) {
  result <- fastcpd(
    data = data.frame(y = data[, 1], x = data[, -1]), family = "lasso", ...
  )
  result@call <- match.call()
  result
}

#' @rdname fastcpd_lasso
#' @export
fastcpd.lasso <- fastcpd_lasso  # nolint: Conventional R function style

#' @title Find change points efficiently in linear regression models
#' @param data A matrix or a data frame with the response variable as the first
#' column.
#' @param ... Other arguments passed to [fastcpd()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [fastcpd_lm()] and [fastcpd.lm()] are wrapper
#' functions of [fastcpd()] to find change points in linear
#' regression models. The function is similar to [fastcpd()] except that
#' the data is by default a matrix or data frame with the response variable
#' as the first column and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_lm_1.R
#' @example tests/testthat/examples/fastcpd_lm_2.txt
#' @seealso [fastcpd()]
#'
#' @md
#' @rdname fastcpd_lm
#' @export
fastcpd_lm <- function(data, ...) {
  result <- fastcpd(
    data = data.frame(y = data[, 1], x = data[, -1]), family = "lm", ...
  )
  result@call <- match.call()
  result
}

#' @rdname fastcpd_lm
#' @export
fastcpd.lm <- fastcpd_lm  # nolint: Conventional R function style

#' @title Find change points efficiently in mean change models
#' @param data A matrix, a data frame or a vector.
#' @param ... Other arguments passed to [fastcpd()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [fastcpd_mean()] and [fastcpd.mean()] are wrapper
#' functions of [fastcpd()] to find the mean change. The function is
#' similar to [fastcpd()] except that the data is by default a matrix or
#' data frame or a vector with each row / element as an observation and thus a
#' formula is not required here.
#' @example tests/testthat/examples/fastcpd_mean_1.R
#' @example tests/testthat/examples/fastcpd_mean_2.R
#' @seealso [fastcpd()]
#'
#' @md
#' @rdname fastcpd_mean
#' @export
fastcpd_mean <- function(data, ...) {
  if (is.null(dim(data)) || length(dim(data)) == 1) {
    data <- matrix(data, ncol = 1)
  }
  result <- fastcpd(
    formula = ~ . - 1, data = data.frame(x = data), family = "mean", ...
  )
  result@call <- match.call()
  result
}

#' @rdname fastcpd_mean
#' @export
fastcpd.mean <- fastcpd_mean  # nolint: Conventional R function style

#' @title Find change points efficiently in mean variance change models
#' @param data A matrix, a data frame or a vector.
#' @param ... Other arguments passed to [fastcpd()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [fastcpd_meanvariance()], [fastcpd.meanvariance()],
#' [fastcpd_mv()], [fastcpd.mv()] are wrapper
#' functions of [fastcpd()] to find the meanvariance change. The
#' function is similar to [fastcpd()] except that the data is by
#' default a matrix or data frame or a vector with each row / element as an
#' observation and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_meanvariance_1.R
#' @example tests/testthat/examples/fastcpd_meanvariance_2.R
#' @seealso [fastcpd()]
#'
#' @md
#' @rdname fastcpd_meanvariance
#' @export
fastcpd_meanvariance <- function(data, ...) {
  if (is.null(dim(data)) || length(dim(data)) == 1) {
    data <- matrix(data, ncol = 1)
  }
  result <- fastcpd(
    formula = ~ . - 1, data = data.frame(x = data), family = "meanvariance", ...
  )
  result@call <- match.call()
  result
}

#' @rdname fastcpd_meanvariance
#' @export
fastcpd.meanvariance <-  # nolint: Conventional R function style
  fastcpd_meanvariance

#' @rdname fastcpd_meanvariance
#' @export
fastcpd_mv <- fastcpd_meanvariance

#' @rdname fastcpd_meanvariance
#' @export
fastcpd.mv <- fastcpd_meanvariance  # nolint: Conventional R function style

#' @title Find change points efficiently in Poisson regression models
#' @param data A matrix or a data frame with the response variable as the first
#' column.
#' @param ... Other arguments passed to [fastcpd()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [fastcpd_poisson()] and [fastcpd.poisson()] are
#' wrapper functions of [fastcpd()] to find change points in
#' Poisson regression models. The function is similar to [fastcpd()]
#' except that the data is by default a matrix or data frame with the response
#' variable as the first column and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_poisson.txt
#' @seealso [fastcpd()]
#'
#' @md
#' @rdname fastcpd_poisson
#' @export
fastcpd_poisson <- function(data, ...) {
  result <- fastcpd(
    data = data.frame(y = data[, 1], x = data[, -1]), family = "poisson", ...
  )
  result@call <- match.call()
  result
}

#' @rdname fastcpd_poisson
#' @export
fastcpd.poisson <- fastcpd_poisson  # nolint: Conventional R function style

#' @title Find change points efficiently in time series data
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param family A character string specifying the family of the time series.
#' The value should be one of \code{"ar"}, \code{"var"}, \code{"arima"} or
#' \code{"garch"}.
#' @param order A positive integer or a vector of length less than four
#' specifying the order of the time series. Possible combinations with
#' \code{family} are:
#' \itemize{
#' \item \code{"ar"}, NUMERIC(1): AR(\eqn{p}) model using linear regression.
#' \item \code{"var"}, NUMERIC(1): VAR(\eqn{p}) model using linear regression.
#' \item \code{"arima"}, NUMERIC(3): ARIMA(\eqn{p}, \eqn{d}, \eqn{q}) model
#'   using [stats::arima()].
#' \item \code{"garch"}, NUMERIC(2): GARCH(\eqn{p}, \eqn{q}) model.
#' }
#' @param ... Other arguments passed to [fastcpd()], for example,
#' \code{segment_count}. One special argument can be passed here is
#' \code{include.mean}, which is a logical value indicating whether the
#' mean should be included in the model. The default value is \code{TRUE}.
#' @return A [fastcpd-class] object.
#' @description [fastcpd_ts()] and [fastcpd.ts()] are wrapper functions for
#' [fastcpd()] to find change points in time series data. The function is
#' similar to [fastcpd()] except that the data is a time series and the
#' family is one of \code{"ar"}, \code{"var"}, \code{"arma"}, \code{"arima"} or
#' \code{"garch"}.
#' @example tests/testthat/examples/fastcpd_ts.txt
#' @seealso [fastcpd()]
#'
#' @md
#' @rdname fastcpd_ts
#' @export
fastcpd_ts <- function(data, family = NULL, order = c(0, 0, 0), ...) {
  if (!is.null(family)) {
    family <- tolower(family)
  }

  check_family(family, c("ar", "var", "arima", "arma", "garch"))
  stopifnot(check_order(order, family))

  # TODO(doccstat): Deal with different data types.
  result <- fastcpd(
    formula = ~ . - 1,
    data = data.frame(x = data),
    family = family,
    order = order,
    ...
  )
  result@call <- match.call()
  result
}

#' @rdname fastcpd_ts
#' @export
fastcpd.ts <- fastcpd_ts  # nolint: Conventional R function style

#' @title Find change points efficiently in VAR(\eqn{p}) models
#' @param data A matrix, a data frame or a time series object.
#' @param order A positive integer specifying the order of the VAR model.
#' @param ... Other arguments passed to [fastcpd()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [fastcpd_var()] and [fastcpd.var()] are
#' wrapper functions of [fastcpd_ts()] to find change points in
#' VAR(\eqn{p}) models. The function is similar to [fastcpd_ts()]
#' except that the data is by default a matrix with row as an observation
#' and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_var.R
#' @seealso [fastcpd()]
#'
#' @md
#' @rdname fastcpd_var
#' @export
fastcpd_var <- function(data, order = 0, ...) {
  result <- fastcpd.ts(data, "var", order, ...)
  result@call <- match.call()
  result
}

#' @rdname fastcpd_var
#' @export
fastcpd.var <- fastcpd_var  # nolint: Conventional R function style

#' @title Find change points efficiently in variance change models
#' @param data A matrix, a data frame or a vector.
#' @param ... Other arguments passed to [fastcpd()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [fastcpd_variance()] and [fastcpd.variance()] are wrapper
#' functions of [fastcpd()] to find the variance change. The
#' function is similar to [fastcpd()] except that the data is by
#' default a matrix or data frame or a vector with each row / element as an
#' observation and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_variance_1.R
#' @example tests/testthat/examples/fastcpd_variance_2.R
#' @seealso [fastcpd()]
#'
#' @md
#' @rdname fastcpd_variance
#' @export
fastcpd_variance <- function(data, ...) {
  if (is.null(dim(data)) || length(dim(data)) == 1) {
    data <- matrix(data, ncol = 1)
  }
  result <- fastcpd(
    formula = ~ . - 1, data = data.frame(x = data), family = "variance", ...
  )
  result@call <- match.call()
  result
}

#' @rdname fastcpd_variance
#' @export
fastcpd.variance <- fastcpd_variance  # nolint: Conventional R function style
