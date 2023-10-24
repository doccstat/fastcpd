#' @title Find change points efficiently in time series data
#'
#' @description `fastcpd_ts` is a wrapper function for `fastcpd` to find
#'   change points in time series data. The function is similar to `fastcpd`
#'   except that the data is a time series data and the family is one of
#'   \code{"ar"}, \code{"var"}, \code{"arima"} or \code{"garch"}.
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
