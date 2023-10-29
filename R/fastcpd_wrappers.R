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
fastcpd.binomial <- function(data, ...) {  # nolint: Conventional R function style
  fastcpd(
    data = data.frame(y = data[, 1], x = data[, -1]), family = "binomial", ...
  )
}

#' @rdname fastcpd_binomial
#' @export
fastcpd_binomial <- fastcpd.binomial

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
#'   \code{segment_count}. One special argument can be passed here is
#'   \code{include.mean}, which is a logical value indicating whether the
#'   mean should be included in the model. The default value is \code{TRUE}.
#'
#' @return A class \code{fastcpd} object.
#'
#' @rdname fastcpd_garch
#' @export
fastcpd.garch <- function(  # nolint: Conventional R function style
  data,
  order = c(0, 0),
  ...
) {
  fastcpd.ts(c(data), "garch", order, ...)
}

#' @rdname fastcpd_garch
#' @export
fastcpd_garch <- fastcpd.garch

#' @title Find change points efficiently in penalized linear regression models
#'
#' @description \code{"fastcpd_lasso"} and \code{"fastcpd.lasso"} are wrapper
#'   functions of \code{\link{fastcpd}} to find change points in penalized
#'   linear regression models. The function is similar to \code{"fastcpd"}
#'   except that the data is by default a matrix or data frame with the response
#'   variable as the first column and thus a formula is not required here.
#'
#' @example tests/testthat/examples/fastcpd_lasso.R
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
fastcpd.lasso <- function(data, ...) {  # nolint: Conventional R function style
  fastcpd(
    data = data.frame(y = data[, 1], x = data[, -1]), family = "lasso", ...
  )
}

#' @rdname fastcpd_lasso
#' @export
fastcpd_lasso <- fastcpd.lasso

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
fastcpd.lm <- function(data, ...) {  # nolint: Conventional R function style
  fastcpd(data = data.frame(y = data[, 1], x = data[, -1]), family = "lm", ...)
}

#' @rdname fastcpd_lm
#' @export
fastcpd_lm <- fastcpd.lm

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
fastcpd.ma <- function(  # nolint: Conventional R function style
  data,
  order = 0,
  ...
) {
  fastcpd.ts(c(data), "ma", order, ...)
}

#' @rdname fastcpd_ma
#' @export
fastcpd_ma <- fastcpd.ma

#' @title Find change points efficiently in linear regression models
#'
#' @description \code{"fastcpd_mean"} and \code{"fastcpd.mean"} are wrapper
#'   functions of \code{\link{fastcpd}} to find the mean change. The function is
#'   similar to \code{"fastcpd"} except that the data is by default a matrix or
#'   data frame or a vector with each row / element as an observation and thus a
#'   formula is not required here.
#'
#' @example tests/testthat/examples/fastcpd_mean.txt
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
fastcpd.mean <- function(data, ...) {  # nolint: Conventional R function style
  data <- matrix(data)
  p <- ncol(data)
  segment_count <- 10
  if (methods::hasArg("segment_count")) {
    segment_count <- eval.parent(match.call()[["segment_count"]])
  }
  block_size <- max(floor(sqrt(nrow(data)) / (segment_count + 1)), 2)
  block_count <- floor(nrow(data) / block_size)
  data_all_covs <- array(NA, dim = c(block_count, p, p))
  for (block_index in seq_len(block_count)) {
    block_start <- (block_index - 1) * block_size + 1
    block_end <- if (block_index < block_count) {
      block_index * block_size
    } else {
      nrow(data)
    }
    data_all_covs[block_index, , ] <- stats::cov(data[block_start:block_end, ])
  }
  data_all_cov <- colMeans(data_all_covs)
  fastcpd(
    formula = ~ . - 1,
    data = data.frame(x = data),
    cost = function(data) {
      n <- nrow(data)
      demeaned_data <- sweep(data, 2, colMeans(data))
      n / 2 * (
        log(det(data_all_cov)) + p * log(2 * pi) +
          sum(diag(solve(data_all_cov, crossprod(demeaned_data)))) / n
      )
    },
    beta = (p + 1) * log(nrow(data)) / 2,
    ...
  )
}

#' @rdname fastcpd_mean
#' @export
fastcpd_mean <- fastcpd.mean

#' @title Find change points efficiently in Poisson regression models
#'
#' @description \code{"fastcpd_poisson"} and \code{"fastcpd.poisson"} are
#'   wrapper functions of \code{\link{fastcpd}} to find change points in
#'   Poisson regression models. The function is similar to \code{"fastcpd"}
#'   except that the data is by default a matrix or data frame with the response
#'   variable as the first column and thus a formula is not required here.
#'
#' @example tests/testthat/examples/fastcpd_poisson.R
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
fastcpd.poisson <- function(  # nolint: Conventional R function style
  data,
  ...
) {
  fastcpd(
    data = data.frame(y = data[, 1], x = data[, -1]), family = "poisson", ...
  )
}

#' @rdname fastcpd_poisson
#' @export
fastcpd_poisson <- fastcpd.poisson

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
#'   - ar, NUMERIC(3): ARIMA(p, 0, 0) model using `forecast::Arima`, where
#'     \code{p} is the first element of the vector.
#'   - var, NUMERIC(1): VAR(p) model using linear regression.
#'   - ma, NUMERIC(1): MA(q) model using `forecast::Arima`.
#'   - ma, NUMERIC(3): ARIMA(0, 0, q) model using `forecast::Arima`, where
#'     \code{q} is the third element of the vector.
#'   - arima, NUMERIC(3): ARIMA(p, d, q) model using `forecast::Arima`.
#'   - garch, NUMERIC(2): GARCH(p, q) model using `fGarch::garchFit`.
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
fastcpd.ts <- function(  # nolint: Conventional R function style
  data,
  family = NULL,
  order = c(0, 0, 0),
  ...
) {
  if (!is.null(family)) {
    family <- tolower(family)
  }

  # TODO(doccstat): Split into different families.
  stopifnot(check_family(family, c("ar", "var", "ma", "arima", "garch")))
  stopifnot(check_order(order, family))

  # TODO(doccstat): Deal with different data types.
  fastcpd(
    formula = ~ . - 1,
    data = data.frame(x = data),
    family = family,
    order = order,
    ...
  )
}

#' @rdname fastcpd_ts
#' @export
fastcpd_ts <- fastcpd.ts
