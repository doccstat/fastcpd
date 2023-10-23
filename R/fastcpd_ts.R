#' @title Find change points efficiently in time series data
#'
#' @description `fastcpd_ts` is a wrapper function for `fastcpd` to find
#'   change points in time series data. The function is similar to `fastcpd`
#'   except that the data is a time series data and the family is one of
#'   \code{"ar"}, \code{"var"}, \code{"arima"} or \code{"garch"}.
#'
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param family A character string specifying the family of the time series.
#'   The value should be one of \code{"ar"}, \code{"var"}, \code{"arima"} or
#'   \code{"garch"}.
#' @param order A vector of length 3 or a positive integer specifying the order
#'   of the time series. The convention of a vector of length 3 follows the
#'   ARIMA model. If a single positive integer is provided, for "ar" and "var"
#'   family, the model is specified as AR(p) or VAR(p) respectively and for
#'   "arima" model, the model is specified as MA(q). If a vector of length 3 is
#'   provided, the model is specified as ARIMA(p, d, q). AR(p) model can also
#'   be specified as ARIMA(p, 0, 0). "ar" family uses linear regression to fit
#'   the time series while "arima" family uses `forecast::Arima`.
#' @param ... Other arguments passed to \code{\link{fastcpd}}.
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
  allowed_family <- c("ar", "var", "arima", "garch")

  # TODO(doccstat): Verify the documentation of `fastcpd` and `fastcpd_ts`.
  if (is.null(family)) {
    stop(r"[The family should be one of "ar", "var", "arima" or "garch".]")
  }
  family <- tolower(family)

  if (!(family %in% allowed_family)) {
    error_message <- r"[
The family should be one of "ar", "var", "arima" or "garch",
while the provided family is {family}.]"
    stop(gsub("{family}", family, error_message, fixed = TRUE))
  }

  if (all(order == 0)) {
    stop(r"[The order should be specified as a vector of length 3.]")
  }
  if (any(order < 0) || any(order != floor(order))) {
    stop(r"[The order should be non-negative integers.]")
  }

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
