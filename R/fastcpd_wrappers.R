#' @name fastcpd_family
#' @aliases fastcpd.family
#' @title Wrapper functions for fastcpd
#' @description Wrapper functions for fastcpd to find change points in various
#' models.
#' @seealso [fastcpd.mean()], [fastcpd.variance()], [fastcpd.mv()],
#' [fastcpd.meanvariance()] for basic statistics change models;
#' [fastcpd.lm()], [fastcpd.binomial()], [fastcpd.poisson()],
#' [fastcpd.lasso()] for regression coefficients change models;
#' [fastcpd.ar()], [fastcpd.var()], [fastcpd.ma()], [fastcpd.arima()],
#' [fastcpd.arma()], [fastcpd.garch()] for change in time series models.
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
#' @example tests/testthat/examples/fastcpd_ar.R
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
#' @example tests/testthat/examples/fastcpd_lm.R
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

#' @title Find change points efficiently in MA(\eqn{q}) models
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param order A positive integer or a vector of length three with the first
#' two elements being zeros specifying the order of the MA model.
#' @param ... Other arguments passed to [fastcpd()], for example,
#' \code{segment_count}. One special argument can be passed here is
#' \code{include.mean}, which is a logical value indicating whether the
#' mean should be included in the model. The default value is \code{TRUE}.
#' @return A [fastcpd-class] object.
#' @description [fastcpd_ma()] and [fastcpd.ma()] are
#' wrapper functions of [fastcpd()] to find change points in
#' MA(\eqn{q}) models. The function is similar to [fastcpd()]
#' except that the data is by default a one-column matrix or univariate vector
#' and thus a formula is not required here.
#' @example tests/testthat/examples/fastcpd_ma.txt
#' @seealso [fastcpd()]
#'
#' @md
#' @rdname fastcpd_ma
#' @export
fastcpd_ma <- function(data, order = 0, ...) {
  result <- fastcpd.ts(c(data), "ma", order, ...)
  result@call <- match.call()
  result
}

#' @rdname fastcpd_ma
#' @export
fastcpd.ma <- fastcpd_ma  # nolint: Conventional R function style

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
#' @example tests/testthat/examples/fastcpd_mean.R
#' @example tests/testthat/examples/fastcpd_mean-time.R
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
#' @example tests/testthat/examples/fastcpd_meanvariance.R
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
#' \item \code{"ma"}, NUMERIC(1): MA(\eqn{q}) model using [forecast::Arima()].
#' \item \code{"ma"}, NUMERIC(3): ARIMA(0, 0, \eqn{q}) model using
#'   [forecast::Arima()], where \eqn{q} is the third element of the vector.
#' \item \code{"arima"}, NUMERIC(3): ARIMA(\eqn{p}, \eqn{d}, \eqn{q}) model
#'   using [forecast::Arima()].
#' \item \code{"garch"}, NUMERIC(2): GARCH(\eqn{p}, \eqn{q}) model using
#'   [tseries::garch()].
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

  stopifnot(
    check_family(family, c("ar", "var", "ma", "arima", "arma", "garch"))
  )
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
#' @example tests/testthat/examples/fastcpd_variance.R
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
