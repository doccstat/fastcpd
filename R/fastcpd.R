#' @title Find change points efficiently
#' @aliases fastcpd
#' @param formula A formula object specifying the model to be fitted. The
#' (optional) response variable should be on the LHS of the formula, while the
#' covariates should be on the RHS. The naming of variables used in the formula
#' should be consistent with the column names in the data frame provided in
#' \code{data}. The intercept term should be removed from the formula.
#' The response variable is not needed for mean/variance change models and time
#' series models. By default, an intercept column will be added to the data,
#' similar to the [lm()] function.
#' Thus, it is suggested that users should remove the intercept term by
#' appending \code{- 1} to the formula. Note that the [detect_family] functions
#' do not require a formula input.
#' @param data A data frame of dimension \eqn{T \times d}{T * d} containing the
#' data to be segmented (where each row denotes a data point
#' \eqn{z_t \in \mathbb{R}^d}{z_t in R^d} for \eqn{t = 1, \ldots, T}) is
#' required in the main function, while a matrix or a vector input is also
#' accepted in the [detect_family] functions.
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
#' value. [detect()] performs the vanilla PELT algorithm, and
#' \code{cost_gradient} and \code{cost_hessian} should not be provided since no
#' parameter updating is necessary for vanilla PELT.
#' In the second format, the loss function
#' \eqn{\sum_{i = s}^t l(z_i, \theta)}{sum_{i = s}^t l(z_i, \theta)} is
#' provided, which has to be optimized over the parameter \eqn{\theta} to
#' obtain the cost value. A detailed discussion about the custom cost function
#' usage can be found in the references.
#'
#' \strong{Compiled cost functions.} For performance-sensitive use cases,
#' \code{cost} (as well as \code{cost_gradient} and \code{cost_hessian}) may
#' instead be a pre-compiled C++ function passed as an external pointer
#' (\code{externalptr}), avoiding R-call overhead in the hot loop. Build one
#' with \code{Rcpp::XPtr}: write a function matching one of
#' \code{double cost(arma::mat const& data)} (PELT-style, format one above) or
#' \code{double cost(arma::mat const& data, arma::colvec const& theta)}
#' (SeGD-style, format two above), take its address, and wrap it as
#' \code{xptr <- Rcpp::XPtr<FnPtr>(new FnPtr(&your_cost), TRUE,
#' Rcpp::wrap("fastcpd_cost_pelt"))} (or \code{"fastcpd_cost_sen"} for the
#' two-argument form) -- the tag string is required and checked at runtime.
#' Since external pointers carry no \code{formals}, also set
#' \code{attr(xptr, "fastcpd_cost_arity") <- 1L} (or \code{2L} for the
#' two-argument form) so [detect()] can route it like an R closure of the
#' same arity. Pass \code{xptr} as \code{cost} exactly as you would an R
#' closure. A compiled \code{cost} cannot be combined with
#' \code{cost_gradient} / \code{cost_hessian} (those drive an R-level
#' \code{stats::optim} warm start that requires an R closure for \code{cost}).
#' See the custom model vignette for a complete worked example.
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
#' Like \code{cost}, this may also be a compiled function passed as an
#' \code{externalptr} matching
#' \code{arma::colvec cost_gradient(arma::mat const& data, arma::colvec const& theta)},
#' wrapped via \code{Rcpp::XPtr} and tagged \code{"fastcpd_cost_gradient"}.
#' @param cost_hessian Hessian of the custom loss function. The Hessian function
#' takes two inputs, the first being a matrix representing a segment of the
#' data, similar to the format used in the \code{cost} function, and the second
#' being the parameter that needs to be optimized. The gradient function returns
#' the Hessian of the loss function, i.e.,
#' \eqn{\sum_{i = s}^t \nabla^2 l(z_i, \theta)}{sum_{i = s}^t l''(z_i, \theta)}.
#' Like \code{cost}, this may also be a compiled function passed as an
#' \code{externalptr} matching
#' \code{arma::mat cost_hessian(arma::mat const& data, arma::colvec const& theta)},
#' wrapped via \code{Rcpp::XPtr} and tagged \code{"fastcpd_cost_hessian"}.
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
#' \item \code{include.mean} is retained for call compatibility in ARIMA
#' models but must be \code{FALSE}; R and Python share a zero-mean native
#' ARMA likelihood.
#' \item \code{show.progress} is used to control the progress bar. By default
#' no progress bar is shown. Set \code{show.progress = TRUE} to display a
#' tqdm-format progress bar on stderr showing PELT timestep progress.
#' \item \code{p.response} is used to specify the number of response variables.
#' This parameter is especially useful for linear models with multivariate
#' responses.
#' }
#' @return A [fastcpd-class] object.
#' @description [detect()] takes in formulas, data, families and extra
#' parameters and returns a [fastcpd-class] object.
#' @section Gallery:
#' <https://github.com/doccstat/fastcpd-r/tree/main/tests/testthat/examples>
#' @section References:
#' Xingchi Li, Xianyang Zhang (2026). ``fastcpd: Fast Change Point Detection
#' in R.'' _Journal of Statistical Software_, **116**(6), 1--53.
#' \doi{10.18637/jss.v116.i06}.
#'
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
#' @example tests/testthat/examples/fastcpd_custom_xptr.txt
#' @seealso [detect_family] for the family-specific function;
#' [plot.fastcpd()] for plotting the results,
#' [summary.fastcpd()] for summarizing the results.
#'
#' @md
#' @importFrom methods show
#' @importFrom Rcpp evalCpp
#' @export
#' @useDynLib fastcpd, .registration = TRUE
detect <- function(  # nolint: cyclomatic complexity
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
  trim = 0,
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
      "quantile",  # -> "quantile"
      "mean",  # -> "mean"
      "variance",  # -> "variance"
      "meanvariance",  # -> "meanvariance"
      "arma",  # -> "arma"
      "ar",  # -> "gaussian"
      "var",  # -> "mgaussian"
      "arima",  # -> "arima" (segment-local differencing + native ARMA)
      "garch",  # -> "garch"
      "exponential",  # -> "exponential"
      "kcp",          # -> "mean" (random Fourier feature transform)
      "custom"  # -> "custom"
    )
  )

  # Check the validity of the `cost` parameter.
  check_cost(cost, cost_gradient, cost_hessian, family)

  # Check the validity of the `cost_adjustment` parameter.
  cost_adjustment_missing <- missing(cost_adjustment)
  if (is.null(cost_adjustment)) {
    cost_adjustment <- "BIC"
  }
  stopifnot(cost_adjustment %in% c("BIC", "MBIC", "MDL"))

  matrix_input <- FALSE
  if (methods::hasArg(".fastcpd_matrix_input")) {
    matrix_input <- isTRUE(
      eval.parent(match.call()[[".fastcpd_matrix_input"]])
    )
  }
  if (matrix_input) {
    data_ <- as.matrix(data)
    stopifnot("`data` must be numeric." = is.numeric(data_))
    storage.mode(data_) <- "double"
    y <- data_[, 1]
  } else {
    # The following code is adapted from the `lm` function from base R.
    match_formula <- match.call(expand.dots = FALSE)
    matched_formula <- match(c("formula", "data"), names(match_formula), 0L)
    match_formula <- match_formula[c(1L, matched_formula)]
    match_formula$drop.unused.levels <- TRUE
    match_formula[[1L]] <- quote(stats::model.frame)
    match_formula <- eval(match_formula, parent.frame())
    y <- stats::model.response(match_formula, "numeric")
    data_ <- cbind(y, stats::model.matrix(formula, data = data))
  }

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
  } else if (family == "arma") {
    stopifnot("Data should be a univariate time series." = ncol(data_) == 1)
    stopifnot(check_order(order, family))
  }

  # Check the parameters passed in the ellipsis.
  include_mean <- FALSE
  p_response <- get_p_response(family, y, data_)
  r_progress <- FALSE
  if (methods::hasArg("include.mean")) {
    include_mean <- eval.parent(match.call()[["include.mean"]])
  }
  if (methods::hasArg("p.response")) {
    p_response <- eval.parent(match.call()[["p.response"]])
  }
  if (methods::hasArg("show.progress")) {
    r_progress <- eval.parent(match.call()[["show.progress"]])
  }
  native_order <- order
  index_offset <- 0L
  make_ar_design <- function(series, ar_order) {
    if (nrow(series) <= ar_order) {
      stop("AR order must be smaller than the number of rows.")
    }
    response <- series[ar_order + seq_len(nrow(series) - ar_order), ]
    predictors <- matrix(NA, nrow(series) - ar_order, ar_order)
    for (order_i in seq_len(ar_order)) {
      predictors[, order_i] <-
        series[(ar_order - order_i) + seq_len(nrow(series) - ar_order), ]
    }
    cbind(response, predictors)
  }

  if (family %in% c("binomial", "poisson", "lasso")) {
    fastcpd_family <- family
  } else if (family == "quantile") {
    fastcpd_family <- "quantile"
    stopifnot(
      "order[1] must be a numeric value in (0, 1) (the quantile level tau)." =
        is.numeric(order) && order[1] > 0 && order[1] < 1
    )
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
  } else if (family == "exponential") {
    fastcpd_family <- family
    vanilla_percentage <- 1
    p <- ncol(data_)
  } else if (family == "kcp") {
    original_p <- ncol(data_)
    n <- nrow(data_)
    D <- as.integer(if (length(order) >= 1 && order[1] > 0) order[1] else 100L)
    sigma_kcp <- if (length(order) >= 2) order[2] else 0
    if (sigma_kcp <= 0) {
      idx <- if (n > 1000L) sample.int(n, 1000L) else seq_len(n)
      d2 <- as.numeric(stats::dist(data_[idx, , drop = FALSE])^2)
      sigma_kcp <- sqrt(stats::median(d2[d2 > 0]) / 2)
    }
    omega <- matrix(stats::rnorm(original_p * D, sd = 1 / sigma_kcp), original_p, D)
    b <- stats::runif(D, 0, 2 * pi)
    data_ <- sqrt(2 / D) * cos(data_ %*% omega + matrix(b, n, D, byrow = TRUE))
    fastcpd_family <- "mean"
    vanilla_percentage <- 1
    p <- D
    if (is.character(beta)) beta <- (original_p + 2) * log(n) / 2
    if (is.null(variance_estimation)) variance_estimation <- diag(D)
    if (cost_adjustment_missing) cost_adjustment <- "BIC"
  } else if (family == "garch") {
    p <- sum(order) + 1
    fastcpd_family <- family
    vanilla_percentage <- 1
  } else if (family == "lm" && p_response == 1) {
    fastcpd_family <- "gaussian"
  } else if (family == "ar") {
    p <- order
    fastcpd_family <- "gaussian"
    data_ <- make_ar_design(data_, p)
    index_offset <- p
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
    index_offset <- order
  } else if (family == "arma" && order[1] == 0) {
    p <- sum(order) + 1
    fastcpd_family <- "ma"
  } else if (family == "arma" && order[2] == 0) {
    p <- order[1]
    fastcpd_family <- "gaussian"
    data_ <- make_ar_design(data_, p)
    index_offset <- p
  } else if (family == "arma" && order[1] != 0) {
    p <- sum(order) + 1
    fastcpd_family <- family
  } else if (family == "arima") {
    stopifnot(
      "`include.mean` must be a single non-missing logical value." =
        is.logical(include_mean) && length(include_mean) == 1 &&
          !is.na(include_mean)
    )
    if (include_mean) {
      stop(
        "`include.mean = TRUE` is not supported by the unified ARIMA ",
        "likelihood; use the default `FALSE`."
      )
    }
    p <- sum(order[-2]) + 1
    if (order[2] == 0) {
      native_order <- order[c(1, 3)]
      if (order[1] == 0) {
        fastcpd_family <- "ma"
      } else if (order[3] == 0) {
        p <- order[1]
        fastcpd_family <- "gaussian"
        data_ <- make_ar_design(data_, p)
        index_offset <- p
      } else {
        fastcpd_family <- "arma"
      }
    } else {
      fastcpd_family <- "arima"
      vanilla_percentage <- 1
    }
  } else {
    if (!methods::hasArg("p")) {
      p <- ncol(data_) - 1
    }
    fastcpd_family <- "custom"
    if (!is.null(cost) && cost_arity(cost) == 1) {
      vanilla_percentage <- 1
    }
  }

  cost_pelt <- NULL
  cost_sen <- NULL
  if (!is.null(cost) && cost_arity(cost) == 1) {
    cost_pelt <- cost
  } else {
    cost_sen <- cost
  }

  sigma_ <- if (!is.null(variance_estimation)) {
    as.matrix(variance_estimation)
  } else if (family == "mean") {
    variance.mean(data_)
  } else if (family == "var" || family == "lm" && p_response > 1) {
    nearest_pd_(variance.lm(data_, p_response))
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

  # Disable pruning for families whose recursive/non-additive costs do not
  # satisfy PELT's pruning condition. Adjust the pruning coefficient for MBIC
  # and MDL.
  pruning_coef <- get_pruning_coef(
    methods::hasArg("pruning_coef"),
    pruning_coef,
    cost_adjustment,
    fastcpd_family,
    p
  )

  result <- fastcpd_impl(
    data_, beta, cost_adjustment, segment_count, trim, momentum_coef,
    multiple_epochs, fastcpd_family, epsilon, p, native_order, cost_pelt,
    cost_sen,
    cost_gradient, cost_hessian, cp_only, vanilla_percentage, warm_start,
    lower, upper, line_search, sigma_, p_response, pruning_coef, r_progress
  )

  raw_cp_set <- c(result$raw_cp_set)
  cp_set <- c(result$cp_set)

  if (index_offset > 0) {
    raw_cp_set <- raw_cp_set + index_offset
    cp_set <- cp_set + index_offset
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
      expr = if (index_offset > 0 && family != "var") {
        residuals <- matrix(c(rep(NA, index_offset), residuals))
      } else if (family == "var") {
        residuals <- rbind(
          matrix(NA, nrow = order, ncol = ncol(residuals)), residuals
        )
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

#' @rdname detect
#' @export
fastcpd <- detect

#' @name detect_family
#' @aliases fastcpd_family fastcpd.family
#' @title Family-specific change point detection
#' @description Convenience functions built on [detect()] for common models.
#' @seealso
#' Basic statistics:
#' [detect_mean()], [detect_variance()], [detect_meanvariance()],
#' [detect_exponential()];
#'
#' Regression:
#' [detect_lm()], [detect_binomial()], [detect_poisson()],
#' [detect_lasso()], [detect_quantile()];
#'
#' Time series:
#' [detect_ar()], [detect_var()], [detect_arima()], [detect_arma()],
#' [detect_garch()];
#'
#' Distribution-free:
#' [detect_kernel()], [detect_rank()].
#'
#' @md
NULL

#' @title Find change points efficiently in AR(\eqn{p}) models
#' @aliases fastcpd_ar fastcpd.ar
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param order A positive integer specifying the order of the AR model.
#' @param ... Other arguments passed to [detect()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [detect_ar()] and [fastcpd.ar()] are
#' wrapper functions of [detect()] to find change points in
#' AR(\eqn{p}) models. The function is similar to [detect()] except that
#' the data is by default a one-column matrix or univariate vector
#' and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_ar_1.R
#' @example tests/testthat/examples/fastcpd_ar_2.R
#' @seealso [detect()]
#'
#' @md
#' @rdname detect_ar
#' @export
detect_ar <- function(data, order = 0, ...) {
  result <- detect(
    formula = ~ . - 1,
    data = data.frame(x = c(data)),
    family = "ar",
    order = order,
    ...
  )
  result@call <- match.call()
  result
}

#' @rdname detect_ar
#' @export
fastcpd_ar <- detect_ar

#' @rdname detect_ar
#' @export
fastcpd.ar <- detect_ar  # nolint: Conventional R function style

#' @title Find change points efficiently in
#' ARIMA(\eqn{p}, \eqn{d}, \eqn{q}) models
#' @aliases fastcpd_arima fastcpd.arima
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param order A vector of length three specifying the order of the ARIMA
#' model.
#' @param include.mean Must be \code{FALSE}. The unified R/Python likelihood
#'   is a zero-mean ARMA likelihood applied after differencing each candidate
#'   segment independently.
#' @param ... Other arguments passed to [detect()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [detect_arima()] and [fastcpd.arima()] are
#' wrapper functions of [detect()] to find change points in
#' ARIMA(\eqn{p}, \eqn{d}, \eqn{q}) models.
#' Differencing is performed independently inside each candidate segment, so
#' no artificial difference is formed across a proposed change-point boundary.
#' Change-point indices refer to the original, undifferenced series. When
#' \eqn{d = 0}, this is the same model as [detect_arma()].
#' The function is similar to [detect()]
#' except that the data is by default a one-column matrix or univariate vector
#' and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_arima.txt
#' @seealso [detect()]
#'
#' @md
#' @rdname detect_arima
#' @export
detect_arima <- function(
  data,
  order = c(1, 1, 0),
  include.mean = FALSE,
  ...
) {
  result <- detect(
    formula = ~ . - 1,
    data = data.frame(x = c(data)),
    family = "arima",
    order = order,
    include.mean = include.mean,
    ...
  )
  result@call <- match.call()
  result
}

#' @rdname detect_arima
#' @export
fastcpd_arima <- detect_arima

#' @rdname detect_arima
#' @export
fastcpd.arima <- detect_arima  # nolint: Conventional R function style

#' @title Find change points efficiently in ARMA(\eqn{p}, \eqn{q}) models
#' @aliases fastcpd_arma fastcpd.arma
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param order A vector of length two specifying the order of the ARMA
#' model.
#' @param ... Other arguments passed to [detect()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [detect_arma()] and [fastcpd.arma()] are
#' wrapper functions of [detect()] to find change points in
#' ARMA(\eqn{p}, \eqn{q}) models. The function is similar to [detect()]
#' except that the data is by default a one-column matrix or univariate vector
#' and thus a formula is not required here. Pure AR models (\eqn{q = 0}) use
#' the same lagged Gaussian path as [detect_ar()]; pure MA models
#' (\eqn{p = 0}) use the native MA family.
#' @example tests/testthat/examples/fastcpd_arma.txt
#' @seealso [detect()]
#'
#' @md
#' @rdname detect_arma
#' @export
detect_arma <- function(data, order = c(0, 0), ...) {
  result <- detect(
    formula = ~ . - 1,
    data = data.frame(x = c(data)),
    family = "arma",
    order = order,
    ...
  )
  result@call <- match.call()
  result
}

#' @rdname detect_arma
#' @export
fastcpd_arma <- detect_arma

#' @rdname detect_arma
#' @export
fastcpd.arma <- detect_arma  # nolint: Conventional R function style

#' @title Find change points efficiently in logistic regression models
#' @aliases detect_logistic_regression fastcpd_binomial fastcpd.binomial
#' @param data A matrix or a data frame with the response variable as the first
#' column.
#' @param ... Other arguments passed to [detect()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [detect_binomial()] and [fastcpd.binomial()] are
#' wrapper functions of [detect()] to find change points in
#' logistic regression models. The function is similar to [detect()]
#' except that the data is by default a matrix or data frame with the response
#' variable as the first column and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_binomial.R
#' @seealso [detect()]
#'
#' @md
#' @rdname detect_binomial
#' @export
detect_binomial <- function(data, ...) {
  data_fast <- fastcpd_matrix_fast_path_data(data)
  if (!is.null(data_fast) && !methods::hasArg("formula")) {
    result <- detect(
      data = data_fast, family = "binomial",
      .fastcpd_matrix_input = TRUE, ...
    )
  } else {
    result <- detect(
      data = data.frame(y = data[, 1], x = data[, -1]), family = "binomial", ...
    )
  }
  result@call <- match.call()
  result
}

#' @rdname detect_binomial
#' @export
detect_logistic_regression <- detect_binomial

#' @rdname detect_binomial
#' @export
fastcpd_binomial <- detect_binomial

#' @rdname detect_binomial
#' @export
fastcpd.binomial <- detect_binomial  # nolint: Conventional R function style

#' @title Find change points efficiently in GARCH(\eqn{p}, \eqn{q}) models
#' @aliases fastcpd_garch fastcpd.garch
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param order A positive integer vector of length two specifying the order of
#' the GARCH model.
#' @param ... Other arguments passed to [detect()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [detect_garch()] and [fastcpd.garch()] are
#' wrapper functions of [detect()] to find change points in
#' GARCH(\eqn{p}, \eqn{q}) models. The function is similar to [detect()]
#' except that the data is by default a one-column matrix or univariate vector
#' and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_garch.txt
#' @example tests/testthat/examples/fastcpd_garch_2.txt
#' @example tests/testthat/examples/fastcpd_garch_3.txt
#' @seealso [detect()]
#'
#' @md
#' @rdname detect_garch
#' @export
detect_garch <- function(data, order = c(0, 0), ...) {
  result <- detect(
    formula = ~ . - 1,
    data = data.frame(x = c(data)),
    family = "garch",
    order = order,
    ...
  )
  result@call <- match.call()
  result
}

#' @rdname detect_garch
#' @export
fastcpd_garch <- detect_garch

#' @rdname detect_garch
#' @export
fastcpd.garch <- detect_garch  # nolint: Conventional R function style

#' @title Find change points efficiently in penalized linear regression models
#' @aliases fastcpd_lasso fastcpd.lasso
#' @param data A matrix or a data frame with the response variable as the first
#' column.
#' @param ... Other arguments passed to [detect()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [detect_lasso()] and [fastcpd.lasso()] are wrapper
#' functions of [detect()] to find change points in penalized
#' linear regression models. The function is similar to [detect()]
#' except that the data is by default a matrix or data frame with the response
#' variable as the first column and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_lasso.txt
#' @seealso [detect()]
#'
#' @md
#' @rdname detect_lasso
#' @export
detect_lasso <- function(data, ...) {
  data_fast <- fastcpd_matrix_fast_path_data(data)
  if (!is.null(data_fast) && !methods::hasArg("formula")) {
    result <- detect(
      data = data_fast, family = "lasso",
      .fastcpd_matrix_input = TRUE, ...
    )
  } else {
    result <- detect(
      data = data.frame(y = data[, 1], x = data[, -1]), family = "lasso", ...
    )
  }
  result@call <- match.call()
  result
}

#' @rdname detect_lasso
#' @export
fastcpd_lasso <- detect_lasso

#' @rdname detect_lasso
#' @export
fastcpd.lasso <- detect_lasso  # nolint: Conventional R function style

#' @title Find change points efficiently in linear regression models
#' @aliases detect_linear_regression fastcpd_lm fastcpd.lm
#' @param data A matrix or a data frame with the response variable as the first
#' column.
#' @param ... Other arguments passed to [detect()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [detect_lm()] and [fastcpd.lm()] are wrapper
#' functions of [detect()] to find change points in linear
#' regression models. The function is similar to [detect()] except that
#' the data is by default a matrix or data frame with the response variable
#' as the first column and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_lm_1.R
#' @example tests/testthat/examples/fastcpd_lm_2.txt
#' @seealso [detect()]
#'
#' @md
#' @rdname detect_lm
#' @export
detect_lm <- function(data, ...) {
  data_fast <- fastcpd_matrix_fast_path_data(data)
  if (!is.null(data_fast) && !methods::hasArg("formula")) {
    result <- detect(
      data = data_fast, family = "lm",
      .fastcpd_matrix_input = TRUE, ...
    )
  } else {
    result <- detect(
      data = data.frame(y = data[, 1], x = data[, -1]), family = "lm", ...
    )
  }
  result@call <- match.call()
  result
}

#' @rdname detect_lm
#' @export
detect_linear_regression <- detect_lm

#' @rdname detect_lm
#' @export
fastcpd_lm <- detect_lm

#' @rdname detect_lm
#' @export
fastcpd.lm <- detect_lm  # nolint: Conventional R function style

#' @title Find change points efficiently in mean change models
#' @aliases fastcpd_mean fastcpd.mean
#' @param data A matrix, a data frame or a vector.
#' @param ... Other arguments passed to [detect()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [detect_mean()] and [fastcpd.mean()] are wrapper
#' functions of [detect()] to find the mean change. The function is
#' similar to [detect()] except that the data is by default a matrix or
#' data frame or a vector with each row / element as an observation and thus a
#' formula is not required here.
#' @example tests/testthat/examples/fastcpd_mean_1.R
#' @example tests/testthat/examples/fastcpd_mean_2.R
#' @seealso [detect()]
#'
#' @md
#' @rdname detect_mean
#' @export
detect_mean <- function(data, ...) {
  data_fast <- fastcpd_matrix_fast_path_data(data, vector_ok = TRUE)
  if (!is.null(data_fast)) {
    result <- detect(
      formula = ~ . - 1, data = data_fast, family = "mean",
      .fastcpd_matrix_input = TRUE, ...
    )
  } else {
    result <- detect(
      formula = ~ . - 1, data = data.frame(x = data), family = "mean", ...
    )
  }
  result@call <- match.call()
  result
}

#' @rdname detect_mean
#' @export
fastcpd_mean <- detect_mean

#' @rdname detect_mean
#' @export
fastcpd.mean <- detect_mean  # nolint: Conventional R function style

#' @title Find change points efficiently in exponentially distributed data
#' @aliases fastcpd_exponential fastcpd.exponential
#' @param data A matrix, a data frame or a vector.
#' @param ... Other arguments passed to [detect()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [detect_exponential()] and [fastcpd.exponential()] are
#' wrapper functions of [detect()] to find changes in the rate of
#' exponentially distributed data, i.e. mean change under exponentially
#' distributed noise (cf. \code{changepoint::cpt.meanvar} with
#' \code{test.stat = "Exponential"}). The function is similar to
#' [detect()] except that the data is by default a matrix or data frame or
#' a vector with each row / element as an observation and thus a formula is
#' not required here.
#' @example tests/testthat/examples/fastcpd_exponential_1.R
#' @seealso [detect()]
#'
#' @md
#' @rdname detect_exponential
#' @export
detect_exponential <- function(data, ...) {
  data_fast <- fastcpd_matrix_fast_path_data(data, vector_ok = TRUE)
  if (!is.null(data_fast)) {
    result <- detect(
      formula = ~ . - 1, data = data_fast, family = "exponential",
      .fastcpd_matrix_input = TRUE, ...
    )
  } else {
    result <- detect(
      formula = ~ . - 1, data = data.frame(x = data), family = "exponential",
      ...
    )
  }
  result@call <- match.call()
  result
}

#' @rdname detect_exponential
#' @export
fastcpd_exponential <- detect_exponential

#' @rdname detect_exponential
#' @export
fastcpd.exponential <-  # nolint: Conventional R function style
  detect_exponential

#' @title Find change points efficiently in mean variance change models
#' @aliases detect_mean_variance fastcpd_meanvariance fastcpd.meanvariance fastcpd_mv fastcpd.mv
#' @param data A matrix, a data frame or a vector.
#' @param ... Other arguments passed to [detect()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [detect_meanvariance()], [fastcpd.meanvariance()],
#' [fastcpd_mv()], [fastcpd.mv()] are wrapper
#' functions of [detect()] to find the meanvariance change. The
#' function is similar to [detect()] except that the data is by
#' default a matrix or data frame or a vector with each row / element as an
#' observation and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_meanvariance_1.R
#' @example tests/testthat/examples/fastcpd_meanvariance_2.R
#' @seealso [detect()]
#'
#' @md
#' @rdname detect_meanvariance
#' @export
detect_meanvariance <- function(data, ...) {
  data_fast <- fastcpd_matrix_fast_path_data(data, vector_ok = TRUE)
  if (!is.null(data_fast)) {
    result <- detect(
      formula = ~ . - 1, data = data_fast, family = "meanvariance",
      .fastcpd_matrix_input = TRUE, ...
    )
  } else {
    result <- detect(
      formula = ~ . - 1, data = data.frame(x = data), family = "meanvariance",
      ...
    )
  }
  result@call <- match.call()
  result
}

#' @rdname detect_meanvariance
#' @export
detect_mean_variance <- detect_meanvariance

#' @rdname detect_meanvariance
#' @export
fastcpd_meanvariance <- detect_meanvariance

#' @rdname detect_meanvariance
#' @export
fastcpd.meanvariance <-  # nolint: Conventional R function style
  detect_meanvariance

#' @rdname detect_meanvariance
#' @export
fastcpd_mv <- detect_meanvariance

#' @rdname detect_meanvariance
#' @export
fastcpd.mv <- detect_meanvariance  # nolint: Conventional R function style

#' @title Find change points efficiently in Poisson regression models
#' @aliases detect_poisson_regression fastcpd_poisson fastcpd.poisson
#' @param data A matrix or a data frame with the response variable as the first
#' column.
#' @param ... Other arguments passed to [detect()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [detect_poisson()] and [fastcpd.poisson()] are
#' wrapper functions of [detect()] to find change points in
#' Poisson regression models. The function is similar to [detect()]
#' except that the data is by default a matrix or data frame with the response
#' variable as the first column and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_poisson.txt
#' @seealso [detect()]
#'
#' @md
#' @rdname detect_poisson
#' @export
detect_poisson <- function(data, ...) {
  data_fast <- fastcpd_matrix_fast_path_data(data)
  if (!is.null(data_fast) && !methods::hasArg("formula")) {
    result <- detect(
      data = data_fast, family = "poisson",
      .fastcpd_matrix_input = TRUE, ...
    )
  } else {
    result <- detect(
      data = data.frame(y = data[, 1], x = data[, -1]), family = "poisson", ...
    )
  }
  result@call <- match.call()
  result
}

#' @rdname detect_poisson
#' @export
detect_poisson_regression <- detect_poisson

#' @rdname detect_poisson
#' @export
fastcpd_poisson <- detect_poisson

#' @rdname detect_poisson
#' @export
fastcpd.poisson <- detect_poisson  # nolint: Conventional R function style

#' @title Find change points efficiently in quantile regression models
#' @aliases detect_quantile_regression fastcpd_quantile fastcpd.quantile
#' @param data A matrix or a data frame with the response variable as the first
#' column and covariates in the remaining columns.
#' @param order Quantile level \eqn{\tau}, a numeric value in (0, 1). The
#' default is 0.5, corresponding to median regression.
#' @param ... Other arguments passed to [detect()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [detect_quantile()] and [fastcpd.quantile()] are wrapper
#' functions of [detect()] to detect change points in quantile regression
#' models using the pinball (check function) loss
#' \eqn{\rho_\tau(u) = u(\tau - \mathbf{1}_{u < 0})}.
#' The function detects changes in the conditional \eqn{\tau}-quantile of
#' the response given the covariates. The segment cost is minimised via
#' iteratively reweighted least squares (IRLS).
#' @example tests/testthat/examples/fastcpd_quantile.txt
#' @example tests/testthat/examples/fastcpd_quantile_2.R
#' @seealso [detect()]
#'
#' @md
#' @rdname detect_quantile
#' @export
detect_quantile <- function(data, order = 0.5, ...) {
  stopifnot(
    "order must be a single numeric value in (0, 1) (the quantile level)." =
      is.numeric(order) && length(order) == 1 && order > 0 && order < 1
  )
  data_fast <- fastcpd_matrix_fast_path_data(data)
  if (!is.null(data_fast) && !methods::hasArg("formula")) {
    result <- detect(
      data = data_fast,
      family = "quantile",
      order = order,
      .fastcpd_matrix_input = TRUE,
      ...
    )
  } else {
    result <- detect(
      data = data.frame(y = data[, 1], x = data[, -1]),
      family = "quantile",
      order = order,
      ...
    )
  }
  result@call <- match.call()
  result
}

#' @rdname detect_quantile
#' @export
detect_quantile_regression <- detect_quantile

#' @rdname detect_quantile
#' @export
fastcpd_quantile <- detect_quantile

#' @rdname detect_quantile
#' @export
fastcpd.quantile <- detect_quantile  # nolint: Conventional R function style

#' @title Find change points efficiently via kernel change point detection
#' @aliases detect_kcp fastcpd_kcp fastcpd.kcp
#' @param data A numeric vector or a matrix with one row per observation.
#' @param order A numeric vector of length up to 2: \code{order[1]} is the
#'   number of random Fourier features (default 100; larger values improve
#'   accuracy at the cost of speed) and \code{order[2]} is the RBF kernel
#'   bandwidth \eqn{\sigma} (default 0, which triggers the median heuristic
#'   \eqn{\sigma = \sqrt{\mathrm{median}(d^2)/2}} on a subsample of up to
#'   1000 observations). Call \code{set.seed()} before this function for
#'   reproducible results.
#' @param ... Other arguments passed to [detect()], for example,
#'   \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [detect_kernel()] and [fastcpd.kcp()] detect change points in
#'   the full distribution of the data using the kernel change point (KCP)
#'   approach. The RBF kernel
#'   \eqn{k(x,y)=\exp(-\|x-y\|^2/(2\sigma^2))} is approximated by random
#'   Fourier features via Bochner's theorem, and change points in the
#'   resulting kernel mean embedding are detected with the existing PELT
#'   infrastructure. This is \eqn{O(nD)} in time — compared to
#'   \eqn{O(n^2)} for the exact kernel PELT — where \eqn{D} is
#'   \code{order[1]}.
#' @example tests/testthat/examples/fastcpd_kcp.R
#' @seealso [detect()]
#' @md
#' @rdname detect_kernel
#' @export
detect_kernel <- function(data, order = c(100L, 0), ...) {
  result <- detect(
    formula = ~ . - 1,
    data = data.frame(as.matrix(data)),
    family = "kcp",
    order = order,
    ...
  )
  result@call <- match.call()
  result
}

#' @rdname detect_kernel
#' @export
detect_kcp <- detect_kernel

#' @rdname detect_kernel
#' @export
fastcpd_kcp <- detect_kernel

#' @rdname detect_kernel
#' @export
fastcpd.kcp <- detect_kernel  # nolint: Conventional R function style

#' @title Find change points efficiently via rank-based change point detection
#' @aliases fastcpd_rank fastcpd.rank
#' @param data A numeric vector or a matrix with one row per observation.
#' @param ... Other arguments passed to [detect()], for example,
#'   \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [detect_rank()] and [fastcpd.rank()] detect change points
#'   using a rank-based, distribution-free cost. Each column is replaced by
#'   its global rank centred at zero, and change points in the mean of these
#'   centred ranks are detected with the existing PELT infrastructure. The
#'   result is fully deterministic and requires no bandwidth selection.
#'   The method is most powerful for location shifts; for scale-only or
#'   general distributional changes, [detect_kernel()] is preferable.
#' @example tests/testthat/examples/fastcpd_rank.R
#' @seealso [detect()]
#' @md
#' @rdname detect_rank
#' @export
detect_rank <- function(data, ...) {
  data_mat <- as.matrix(data)
  n <- nrow(data_mat)
  centered_ranks <- apply(data_mat, 2, rank) - (n + 1) / 2
  result <- detect(
    formula = ~ . - 1, data = data.frame(centered_ranks), family = "mean", ...
  )
  result@call <- match.call()
  result
}

#' @rdname detect_rank
#' @export
fastcpd_rank <- detect_rank

#' @rdname detect_rank
#' @export
fastcpd.rank <- detect_rank  # nolint: Conventional R function style

#' @title Find change points efficiently in time series data
#' @aliases fastcpd.ts
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
#'   using segment-local differencing and the shared native ARMA likelihood.
#' \item \code{"garch"}, NUMERIC(2): GARCH(\eqn{p}, \eqn{q}) model.
#' }
#' @param ... Other arguments passed to [detect()], for example,
#' \code{segment_count}. For ARIMA models, \code{include.mean} must remain
#' \code{FALSE} so the shared zero-mean R/Python likelihood is used.
#' @return A [fastcpd-class] object.
#' @description [fastcpd_ts()] and [fastcpd.ts()] are compatibility wrapper
#' functions for [detect()] to find change points in time series data. New code
#' should use a family-specific wrapper such as [detect_ar()] or call [detect()]
#' directly. Like [detect()], it accepts the time-series families
#' \code{"ar"}, \code{"var"}, \code{"arma"}, \code{"arima"}, and
#' \code{"garch"}.
#' @example tests/testthat/examples/fastcpd_ts.txt
#' @seealso [detect()]
#'
#' @md
#' @keywords internal
#' @rdname fastcpd_ts
#' @export
fastcpd_ts <- function(data, family = NULL, order = c(0, 0, 0), ...) {
  if (!is.null(family)) {
    family <- tolower(family)
  }

  check_family(family, c("ar", "var", "arima", "arma", "garch"))
  stopifnot(check_order(order, family))

  # TODO(doccstat): Deal with different data types.
  result <- detect(
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
#' @aliases fastcpd_var fastcpd.var
#' @param data A matrix, a data frame or a time series object.
#' @param order A positive integer specifying the order of the VAR model.
#' @param ... Other arguments passed to [detect()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [detect_var()] and [fastcpd.var()] are
#' wrapper functions of [detect()] to find change points in
#' VAR(\eqn{p}) models. The function is similar to [detect()]
#' except that the data is by default a matrix with row as an observation
#' and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_var.R
#' @seealso [detect()]
#'
#' @md
#' @rdname detect_var
#' @export
detect_var <- function(data, order = 0, ...) {
  result <- detect(
    formula = ~ . - 1,
    data = data.frame(x = data),
    family = "var",
    order = order,
    ...
  )
  result@call <- match.call()
  result
}

#' @rdname detect_var
#' @export
fastcpd_var <- detect_var

#' @rdname detect_var
#' @export
fastcpd.var <- detect_var  # nolint: Conventional R function style

#' @title Find change points efficiently in variance change models
#' @aliases fastcpd_variance fastcpd.variance
#' @param data A matrix, a data frame or a vector.
#' @param ... Other arguments passed to [detect()], for example,
#' \code{segment_count}.
#' @return A [fastcpd-class] object.
#' @description [detect_variance()] and [fastcpd.variance()] are wrapper
#' functions of [detect()] to find the variance change. The
#' function is similar to [detect()] except that the data is by
#' default a matrix or data frame or a vector with each row / element as an
#' observation and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_variance_1.R
#' @example tests/testthat/examples/fastcpd_variance_2.R
#' @seealso [detect()]
#'
#' @md
#' @rdname detect_variance
#' @export
detect_variance <- function(data, ...) {
  data_fast <- fastcpd_matrix_fast_path_data(data, vector_ok = TRUE)
  if (!is.null(data_fast)) {
    result <- detect(
      formula = ~ . - 1, data = data_fast, family = "variance",
      .fastcpd_matrix_input = TRUE, ...
    )
  } else {
    result <- detect(
      formula = ~ . - 1, data = data.frame(x = data), family = "variance", ...
    )
  }
  result@call <- match.call()
  result
}

#' @rdname detect_variance
#' @export
fastcpd_variance <- detect_variance

#' @rdname detect_variance
#' @export
fastcpd.variance <- detect_variance  # nolint: Conventional R function style
