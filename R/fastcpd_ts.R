#' @title Find change points efficiently in time series data
#'
#' @description `fastcpd_ts` is a wrapper function for `fastcpd` to find
#'   change points in time series data. The function is similar to `fastcpd`
#'   except that the data is a time series data and the family is one of
#'   \code{"ar"} or \code{"var"}.
#'
#' @param data A numeric vector, a matrix, a data frame or a time series object.
#' @param family A character string specifying the family of the time series.
#'   The value should be one of \code{"ar"} or \code{"var"}.
#' @param order A positive integer specifying the order of the time series.
#' @param ... Other arguments passed to \code{\link{fastcpd}}.
#'
#' @return A class \code{fastcpd} object.
#'
#' @rdname fastcpd_ts
#' @export
fastcpd.ts <- function(  # nolint: Conventional R function style
  data,
  family = NULL,
  order = NULL,
  ...
) {
  allowed_family <- c("ar", "var")

  stopifnot("The family should be one of \"ar\" or \"var\"" = !is.null(family))
  stopifnot("Order of the time series should be specified" = !is.null(order))
  stopifnot(
    "`order` should be a positive integer" = order > 0 && order == floor(order)
  )

  if (!(family %in% allowed_family)) {
    error_message <- r"[
The family should be one of "ar" or "var"
while the provided family is {family}.]"
    stop(gsub("{family}", family, error_message, fixed = TRUE))
  }

  fastcpd(
    formula = ~ . - 1,
    data = data.frame(x = data),
    family = family,
    p = order,
    ...
  )
}

#' @rdname fastcpd_ts
#' @export
fastcpd_ts <- fastcpd.ts
