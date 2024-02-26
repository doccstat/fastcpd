#' @title Find change points efficiently in AR(p) models
#'
#' @description \code{fastcpd_ar} and \code{fastcpd.ar} are
#'   wrapper functions of \code{\link{fastcpd}} to find change points in
#'   AR(p) models. The function is similar to \code{\link{fastcpd}}
#'   except that the data is by default a one-column matrix or univariate vector
#'   and thus a formula is not required here.
#'
#' @example tests/testthat/examples/fastcpd_ar.R
#'
#' @md
#'
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param order A positive integer specifying the order of the AR model.
#' @param ... Other arguments passed to \code{\link{fastcpd}}, for example,
#'   \code{segment_count}. One special argument can be passed here is
#'   \code{include.mean}, which is a logical value indicating whether the
#'   mean should be included in the model. The default value is \code{TRUE}.
#'
#' @return A class \code{fastcpd} object.
#'
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

#' @title Find change points efficiently in ARIMA(p, d, q) models
#'
#' @description \code{fastcpd_arima} and \code{fastcpd.arima} are
#'   wrapper functions of \code{\link{fastcpd}} to find change points in
#'   ARIMA(p, d, q) models. The function is similar to \code{\link{fastcpd}}
#'   except that the data is by default a one-column matrix or univariate vector
#'   and thus a formula is not required here.
#'
#' @example tests/testthat/examples/fastcpd_arima.txt
#'
#' @md
#'
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param order A vector of length three specifying the order of the ARIMA
#'   model.
#' @param ... Other arguments passed to \code{\link{fastcpd}}, for example,
#'   \code{segment_count}. One special argument can be passed here is
#'   \code{include.mean}, which is a logical value indicating whether the
#'   mean should be included in the model. The default value is \code{TRUE}.
#'
#' @return A class \code{fastcpd} object.
#'
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

#' @title Find change points efficiently in ARMA(p, q) models
#'
#' @description \code{fastcpd_arma} and \code{fastcpd.arma} are
#'   wrapper functions of \code{\link{fastcpd}} to find change points in
#'   ARMA(p, q) models. The function is similar to \code{\link{fastcpd}}
#'   except that the data is by default a one-column matrix or univariate vector
#'   and thus a formula is not required here.
#'
#' @example tests/testthat/examples/fastcpd_arma.txt
#'
#' @md
#'
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param order A vector of length two specifying the order of the ARMA
#'   model.
#' @param ... Other arguments passed to \code{\link{fastcpd}}, for example,
#'   \code{segment_count}.
#'
#' @return A class \code{fastcpd} object.
#'
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
#'
#' @description \code{"fastcpd_binomial"} and \code{"fastcpd.binomial"} are
#'   wrapper functions of \code{\link{fastcpd}} to find change points in
#'   logistic regression models. The function is similar to \code{"fastcpd"}
#'   except that the data is by default a matrix or data frame with the response
#'   variable as the first column and thus a formula is not required here.
#'
#' @example tests/testthat/examples/fastcpd_binomial.R
#'
#' @md
#'
#' @param data A matrix or a data frame with the response variable as the first
#'   column.
#' @param ... Other arguments passed to \code{\link{fastcpd}}, for example,
#'   \code{segment_count}.
#'
#' @return A class \code{fastcpd} object.
#'
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

#' @title Find change points efficiently in GARCH(p, q) models
#'
#' @description \code{"fastcpd_garch"} and \code{"fastcpd.garch"} are
#'   wrapper functions of \code{\link{fastcpd}} to find change points in
#'   GARCH(p, q) models. The function is similar to \code{"fastcpd"}
#'   except that the data is by default a one-column matrix or univariate vector
#'   and thus a formula is not required here.
#'
#' @example tests/testthat/examples/fastcpd_garch.txt
#'
#' @md
#'
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param order A positive integer vector of length two specifying the order of
#'   the GARCH model.
#' @param ... Other arguments passed to \code{\link{fastcpd}}, for example,
#'   \code{segment_count}.
#'
#' @return A class \code{fastcpd} object.
#'
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
#'
#' @description \code{"fastcpd_lasso"} and \code{"fastcpd.lasso"} are wrapper
#'   functions of \code{\link{fastcpd}} to find change points in penalized
#'   linear regression models. The function is similar to \code{"fastcpd"}
#'   except that the data is by default a matrix or data frame with the response
#'   variable as the first column and thus a formula is not required here.
#'
#' @example tests/testthat/examples/fastcpd_lasso.txt
#'
#' @md
#'
#' @param data A matrix or a data frame with the response variable as the first
#'   column.
#' @param ... Other arguments passed to \code{\link{fastcpd}}, for example,
#'   \code{segment_count}.
#'
#' @return A class \code{fastcpd} object.
#'
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
#'
#' @description \code{"fastcpd_lm"} and \code{"fastcpd.lm"} are wrapper
#'   functions of \code{\link{fastcpd}} to find change points in linear
#'   regression models. The function is similar to \code{"fastcpd"} except that
#'   the data is by default a matrix or data frame with the response variable
#'   as the first column and thus a formula is not required here.
#'
#' @example tests/testthat/examples/fastcpd_lm.R
#'
#' @md
#'
#' @param data A matrix or a data frame with the response variable as the first
#'   column.
#' @param ... Other arguments passed to \code{\link{fastcpd}}, for example,
#'   \code{segment_count}.
#'
#' @return A class \code{fastcpd} object.
#'
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

#' @title Find change points efficiently in MA(q) models
#'
#' @description \code{fastcpd_ma} and \code{fastcpd.ma} are
#'   wrapper functions of \code{\link{fastcpd}} to find change points in
#'   MA(q) models. The function is similar to \code{\link{fastcpd}}
#'   except that the data is by default a one-column matrix or univariate vector
#'   and thus a formula is not required here.
#'
#' @example tests/testthat/examples/fastcpd_ma.txt
#'
#' @md
#'
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param order A positive integer or a vector of length three with the first
#'   two elements being zeros specifying the order of the MA model.
#' @param ... Other arguments passed to \code{\link{fastcpd}}, for example,
#'   \code{segment_count}. One special argument can be passed here is
#'   \code{include.mean}, which is a logical value indicating whether the
#'   mean should be included in the model. The default value is \code{TRUE}.
#'
#' @return A class \code{fastcpd} object.
#'
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
#'
#' @description \code{"fastcpd_mean"} and \code{"fastcpd.mean"} are wrapper
#'   functions of \code{\link{fastcpd}} to find the mean change. The function is
#'   similar to \code{"fastcpd"} except that the data is by default a matrix or
#'   data frame or a vector with each row / element as an observation and thus a
#'   formula is not required here.
#'
#' @example tests/testthat/examples/fastcpd_mean.R
#' @example tests/testthat/examples/fastcpd_mean-time.R
#'
#' @md
#'
#' @param data A matrix, a data frame or a vector.
#' @param ... Other arguments passed to \code{\link{fastcpd}}, for example,
#'   \code{segment_count}.
#'
#' @return A class \code{fastcpd} object.
#'
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

#' @title Find change points efficiently in variance change models
#'
#' @description \code{fastcpd_meanvariance}, \code{fastcpd.meanvariance},
#'   \code{fastcpd_mv}, \code{fastcpd.mv} are wrapper
#'   functions of \code{\link{fastcpd}} to find the meanvariance change. The
#'   function is similar to \code{\link{fastcpd}} except that the data is by
#'   default a matrix or data frame or a vector with each row / element as an
#'   observation and thus a formula is not required here.
#'
#' @example tests/testthat/examples/fastcpd_meanvariance.R
#'
#' @md
#'
#' @param data A matrix, a data frame or a vector.
#' @param ... Other arguments passed to \code{\link{fastcpd}}, for example,
#'   \code{segment_count}.
#'
#' @return A class \code{fastcpd} object.
#'
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
#'
#' @description \code{"fastcpd_poisson"} and \code{"fastcpd.poisson"} are
#'   wrapper functions of \code{\link{fastcpd}} to find change points in
#'   Poisson regression models. The function is similar to \code{"fastcpd"}
#'   except that the data is by default a matrix or data frame with the response
#'   variable as the first column and thus a formula is not required here.
#'
#' @example tests/testthat/examples/fastcpd_poisson.txt
#'
#' @md
#'
#' @param data A matrix or a data frame with the response variable as the first
#'   column.
#' @param ... Other arguments passed to \code{\link{fastcpd}}, for example,
#'   \code{segment_count}.
#'
#' @return A class \code{fastcpd} object.
#'
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
#'
#' @description `fastcpd_ts` is a wrapper function for `fastcpd` to find
#'   change points in time series data. The function is similar to `fastcpd`
#'   except that the data is a time series data and the family is one of
#'   \code{"ar"}, \code{"var"}, \code{"arima"} or \code{"garch"}.
#'
#' @example tests/testthat/examples/fastcpd_ts.txt
#'
#' @md
#'
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param family A character string specifying the family of the time series.
#'   The value should be one of \code{"ar"}, \code{"var"}, \code{"arima"} or
#'   \code{"garch"}.
#' @param order A positive integer or a vector of length less than four
#'   specifying the order of the time series. Possible combinations with
#'   `family` are:
#'
#'   - ar, NUMERIC(1): AR(p) model using linear regression.
#'   - var, NUMERIC(1): VAR(p) model using linear regression.
#'   - ma, NUMERIC(1): MA(q) model using `forecast::Arima`.
#'   - ma, NUMERIC(3): ARIMA(0, 0, q) model using `forecast::Arima`, where
#'     \code{q} is the third element of the vector.
#'   - arima, NUMERIC(3): ARIMA(p, d, q) model using `forecast::Arima`.
#'   - garch, NUMERIC(2): GARCH(p, q) model using `tseries::garch`.
#'
#' @param ... Other arguments passed to \code{\link{fastcpd}}, for example,
#'   \code{segment_count}. One special argument can be passed here is
#'   \code{include.mean}, which is a logical value indicating whether the
#'   mean should be included in the model. The default value is \code{TRUE}.
#'
#' @return A class \code{fastcpd} object.
#'
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

#' @title Find change points efficiently in VAR(p) models
#'
#' @description \code{fastcpd_var} and \code{fastcpd.var} are
#'   wrapper functions of \code{\link{fastcpd.ts}} to find change points in
#'   VAR(p) models. The function is similar to \code{\link{fastcpd.ts}}
#'   except that the data is by default a matrix with row as an observation
#'   and thus a formula is not required here.
#'
#' @example tests/testthat/examples/fastcpd_var.R
#'
#' @md
#'
#' @param data A matrix, a data frame or a time series object.
#' @param order A positive integer specifying the order of the VAR model.
#' @param ... Other arguments passed to \code{\link{fastcpd}}, for example,
#'   \code{segment_count}.
#'
#' @return A class \code{fastcpd} object.
#'
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
#'
#' @description \code{fastcpd_variance} and \code{fastcpd.variance} are wrapper
#'   functions of \code{\link{fastcpd}} to find the variance change. The
#'   function is similar to \code{\link{fastcpd}} except that the data is by
#'   default a matrix or data frame or a vector with each row / element as an
#'   observation and thus a formula is not required here.
#'
#' @example tests/testthat/examples/fastcpd_variance.R
#'
#' @md
#'
#' @param data A matrix, a data frame or a vector.
#' @param ... Other arguments passed to \code{\link{fastcpd}}, for example,
#'   \code{segment_count}.
#'
#' @return A class \code{fastcpd} object.
#'
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
