#' An S4 class to store the output created with \link{fastcpd}
#'
#' This S4 class stores the output from \link{fastcpd}. A fastcpd object consist
#' of several slots including the call to \link{fastcpd}, the data used, the
#' family of the model, the change points, the cost values, the residuals, the
#' estimated parameters and a boolean indicating whether the model was fitted
#' with only change points or with change points and parameters, which you can
#' select using \code{@}.
#'
#' @slot call The call to \link{fastcpd}.
#' @slot data The data used.
#' @slot family The family of the model.
#' @slot cp_set The change points.
#' @slot cost_values The cost values for each segment.
#' @slot residuals The residuals for each segment.
#' @slot thetas The estimated parameters for each segment.
#' @slot cp_only A boolean indicating whether the model was fitted with only
#'   change points or with change points and parameters.
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
#' @rdname plot.fastcpd
#' @export
plot.fastcpd <- function(x, ...) {
  y <- x@data[, 1]
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = data.frame(x = seq_len(nrow(x@data)), y = y, label = "response"),
      ggplot2::aes(x = x, y = y)
    ) +
    ggplot2::geom_vline(xintercept = x@cp_set, color = "red")
  if (x@family != "custom" && !x@cp_only) {
    p <- p + ggplot2::geom_point(data = data.frame(x = seq_len(nrow(x@data)),
                                                   y = x@residuals,
                                                   label = "residual"),
                                 ggplot2::aes(x = x, y = y)) +
      ggplot2::facet_wrap(c("label"), ncol = 1, scales = "free_y")
  }
  print(p)
  invisible()
}

#' Plot the data and the change points for a \code{fastcpd} object
#' @param x \code{fastcpd} object.
#' @param y Ignored.
#' @param ... Ignored.
#'
#' @rdname plot.fastcpd
#' @export
setMethod("plot", signature(x = "fastcpd", y = "missing"), plot.fastcpd)

#' @method print fastcpd
#' @rdname print.fastcpd
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
#' @rdname print.fastcpd
#' @export
setMethod("print", "fastcpd", print.fastcpd)

#' @method show fastcpd
#' @rdname show.fastcpd
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
#' @rdname show.fastcpd
#' @export
setMethod("show", "fastcpd", show.fastcpd)

#' @method summary fastcpd
#' @rdname summary.fastcpd
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
      print(object@thetas)
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
#' @rdname summary.fastcpd
#' @export
setMethod("summary", "fastcpd", summary.fastcpd)
