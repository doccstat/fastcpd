#' An S4 class to store the output created with \link{fastcpd}
#'
#' This S4 class stores the output from \link{fastcpd}. A fastcpd object consist
#' of several slots including the call to \link{fastcpd}, the data used, the
#' family of the model, the change points, the cost values, the residuals, the
#' estimated parameters and a boolean indicating whether the model was fitted
#' with only change points or with change points and parameters, which you can
#' select using \code{@}.
#'
#' @slot call The call of the function.
#' @slot data The data passed to the function.
#' @slot family The family of the model.
#' @slot cp_set The set of change points.
#' @slot cost_values The cost function values for each segment.
#' @slot residuals The residuals of the model with change points.
#'   Used only for built-in families.
#' @slot thetas The estimated parameters for each segment. Used only for
#'   built-in families.
#' @slot cp_only A boolean indicating whether `fastcpd` was run to return
#'   only the change points or the change points with the estimated parameters
#'   and cost values for each segment.
#' @export
setClass(
  "fastcpd",
  representation(
    call = "language",
    data = "data.frame",
    family = "character",
    cp_set = "numeric",
    cost_values = "numeric",
    residuals = "numeric",
    thetas = "data.frame",
    cp_only = "logical"
  )
)

#' @method plot fastcpd
#' @rdname plot
#' @export
plot.fastcpd <- function(x, ...) {  # nolint: cyclomatic complexity
  # Plot the built in families only.
  stopifnot(
    "Built-in plot only works for built-in families." = x@family != "custom"
  )
  if (x@family == "mean" && ncol(x@data) > 1) {
    stop("Can not plot mean change points with p > 1.")
  }
  if (!require_namespace("ggplot2")) {
    if (utils_menu() == 1) {
      tryCatch(
        expr = install_packages("ggplot2"),
        error = function(e) {
          stop("ggplot2 could not be installed.")
        }
      )
    } else {
      stop("ggplot2 is not installed. No plot is made.")
    }
  }
  if (require_namespace("ggplot2")) {
    y <- x@data[, 1]
    p <- ggplot2::ggplot() +
      ggplot2::geom_vline(xintercept = x@cp_set, color = "red")

    # Draw lines for time series data and points for other data.
    if (x@family %in% c("ar", "ma", "arma", "arima", "garch")) {
      line_or_point <- ggplot2::geom_line
    } else {
      line_or_point <- ggplot2::geom_point
    }
    p <- p + line_or_point(
      data = data.frame(x = seq_len(nrow(x@data)), y = y, label = "response"),
      ggplot2::aes(x = x, y = y)
    )

    if (x@family != "var" && !x@cp_only) {
      p <- p + ggplot2::geom_point(
        data = data.frame(
          x = seq_len(nrow(x@data)),
          y = x@residuals,
          label = "residual"
        ),
        ggplot2::aes(x = x, y = y),
        na.rm = TRUE
      )
      if (ncol(x@data) == 2 || (x@family == "ar" && nrow(x@thetas) == 1)) {
        xend <- c(x@cp_set, nrow(x@data))
        yend <- as.numeric(x@thetas)
        p <- p + ggplot2::geom_point(
          data = data.frame(
            x = seq_len(nrow(x@data)),
            y = x@data[, ncol(x@data)],
            label = "covariate"
          ),
          ggplot2::aes(x = x, y = y)
        )
        p <- p + ggplot2::geom_segment(
          data = data.frame(
            x = c(1, x@cp_set),
            y = yend,
            xend = xend,
            yend = yend,
            label = "coefficient"
          ),
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
          col = "blue"
        )
      }
      p <- p + ggplot2::facet_wrap("label", nrow = 2, scales = "free_y")
    }
    print(p)
  }
  invisible()
}

#' Plot the data and the change points for a \code{fastcpd} object
#' @param x \code{fastcpd} object.
#' @param y Ignored.
#' @param ... Ignored.
#'
#' @return No return value, called for plotting.
#' @rdname plot
#' @export
setMethod("plot", signature(x = "fastcpd", y = "missing"), plot.fastcpd)

#' @method print fastcpd
#' @rdname print
#' @export
print.fastcpd <- function(x, ...) {
  if (length(x@cp_set)) {
    cat("\nChange points:\n")
    print(x@cp_set)
  } else {
    cat("\nNo change points found\n")
  }
  invisible(x)
}

#' Print the call and the change points for a \code{fastcpd} object
#' @param x \code{fastcpd} object.
#' @param ... Ignored.
#'
#' @return Return a (temporarily) invisible copy of the \code{fastcpd} object.
#'   Called primarily for printing the change points in the model.
#' @rdname print
#' @export
setMethod("print", "fastcpd", print.fastcpd)

#' @method show fastcpd
#' @rdname show
#' @export
show.fastcpd <- function(object) {
  cat(
    "\nA fastcpd object.\n",
    "Available methods to evaluate the object are:\n",
    "plot, print, show, summary\n\n",
    sep = ""
  )
  print(object)
}

#' Show the available methods for a \code{fastcpd} object
#' @param object \code{fastcpd} object.
#'
#' @return No return value, called for showing a list of available methods
#'   for a \code{fastcpd} object.
#' @rdname show
#' @export
setMethod("show", "fastcpd", show.fastcpd)

#' @method summary fastcpd
#' @rdname summary
#' @export
summary.fastcpd <- function(object, ...) {
  cat(
    "\nCall:\n",
    paste(deparse(object@call), sep = "\n", collapse = "\n"), "\n\n",
    sep = ""
  )
  if (length(object@cp_set)) {
    cat("Change points:\n")
    cat(object@cp_set, "\n")
    if (object@family != "custom" && !object@cp_only) {
      cat("\nCost values:\n")
      cat(object@cost_values, "\n")
    }
    if (ncol(object@thetas) > 0) {
      cat("\nParameters:\n")
      if (object@family != "lasso") {
        print(object@thetas)
      } else {
        print(Matrix::Matrix(as.matrix(object@thetas), sparse = TRUE))
      }
    }
  } else {
    cat("No change points found\n")
  }
  invisible(object)
}

#' Show the summary of a \code{fastcpd} object
#' @param object \code{fastcpd} object.
#' @param ... Ignored.
#'
#' @return Return a (temporarily) invisible copy of the \code{fastcpd} object.
#'   Called primarily for printing the summary of the model including the
#'   call, the change points, the cost values and the estimated parameters.
#' @rdname summary
#' @export
setMethod("summary", "fastcpd", summary.fastcpd)
