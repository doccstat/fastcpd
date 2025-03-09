#' @title An S4 class to store the output created with [fastcpd()]
#' @description This S4 class stores the output from [fastcpd()] and
#' [fastcpd.family]. A fastcpd object consist
#' of several slots including the call to [fastcpd()], the data used, the
#' family of the model, the change points, the cost values, the residuals, the
#' estimated parameters and a boolean indicating whether the model was fitted
#' with only change points or with change points and parameters, which you can
#' select using \code{@}.
#' @slot call The call of the function.
#' @slot data The data passed to the function.
#' @slot order The order of the time series model.
#' @slot family The family of the model.
#' @slot cp_set The set of change points.
#' @slot cost_values The cost function values for each segment.
#' @slot residuals The residuals of the model with change points.
#' Used only for built-in families.
#' @slot thetas The estimated parameters for each segment. Used only for
#' built-in families.
#' @slot cp_only A boolean indicating whether [fastcpd()] was run to return
#' only the change points or the change points with the estimated parameters
#' and cost values for each segment.
#' @md
#' @export
setClass(
  "fastcpd",
  representation(
    call = "language",
    data = "data.frame",
    order = "numeric",
    family = "character",
    cp_set = "numeric",
    cost_values = "numeric",
    residuals = "matrix",
    thetas = "data.frame",
    cp_only = "logical"
  )
)

#' @method plot fastcpd
#' @rdname plot
#' @export
plot.fastcpd <- function(  # nolint: cyclomatic complexity
  x,
  color_max_count = Inf,
  data_point_alpha = 0.8,
  data_point_linewidth = 0.5,
  data_point_size = 1,
  legend_position = "none",
  panel_background = ggplot2::element_blank(),
  panel_border = ggplot2::element_rect(fill = NA, colour = "grey20"),
  panel_grid_major = ggplot2::element_line(colour = "grey98"),
  panel_grid_minor = ggplot2::element_line(colour = "grey98"),
  segment_separator_alpha = 0.8,
  segment_separator_color = "grey",
  segment_separator_linetype = "dashed",
  strip_background = ggplot2::element_rect(fill = "grey85", colour = "grey20"),
  xlab = NULL,
  ylab = NULL,
  ...
) {
  if (x@family == "custom") {
    message("Built-in plot should only work for built-in families.")
  }
  if (x@family == "mean" && ncol(x@data) > 1) {
    warning("Can not plot mean change points with p > 1.")
    return()
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 is not installed. No plot is made.")
    return()
  }
  n <- nrow(x@data)
  family <- x@family
  change_points <- sort(c(0, x@cp_set, n))
  color <- rep(seq_along(change_points[-1]), diff(change_points))
  color <- as.factor(color %% color_max_count)

  y <- x@data[, 1]
  p <- ggplot2::ggplot() +
    ggplot2::geom_vline(
      xintercept = x@cp_set,
      color = segment_separator_color,
      linetype = segment_separator_linetype,
      alpha = segment_separator_alpha
    )

  # Draw lines for time series data and points for other data.
  if (family %in% c("ar", "arma", "arima", "garch")) {
    y_label <-
      paste0(toupper(family), "(", paste0(x@order, collapse = ", "), ")")
  } else if (family %in% c("mean", "variance", "meanvariance")) {
    y_label <- "data"
  } else {
    y_label <- "data response"
  }

  data_label_color <- data.frame(
    x = seq_len(n), y = y, label = y_label, color = color
  )
  residual_label_color <- data.frame(
    x = seq_len(n),
    y = x@residuals,
    label = "residual",
    color = color
  )
  covariate_label_color <- data.frame(
    x = seq_len(n),
    y = x@data[, ncol(x@data)],
    label = "covariate",
    color = color
  )
  aesthetic_mapping <- ggplot2::aes(x = x, y = y, color = color)

  if (family %in% c("ar", "arma", "arima", "garch")) {
    p <- p + ggplot2::geom_line(
      data = data_label_color, aesthetic_mapping, alpha = data_point_alpha,
      linewidth = data_point_linewidth
    )
  } else {
    p <- p + ggplot2::geom_point(
      data = data_label_color, aesthetic_mapping, alpha = data_point_alpha,
      size = data_point_size
    )
  }

  if (family != "var" && !x@cp_only) {
    p <- p + ggplot2::geom_point(
      data = residual_label_color,
      aesthetic_mapping,
      na.rm = TRUE,
      alpha = data_point_alpha,
      size = data_point_size
    )
    if (ncol(x@data) == 2 || (family == "ar" && nrow(x@thetas) == 1)) {
      xend <- c(x@cp_set, n)
      yend <- as.numeric(x@thetas)

      coefficient_label <- data.frame(
        x = c(1, x@cp_set),
        y = yend,
        xend = xend,
        yend = yend,
        label = "coefficient"
      )

      p <- p + ggplot2::geom_point(
        data = covariate_label_color,
        aesthetic_mapping,
        size = data_point_size
      )
      p <- p + ggplot2::geom_segment(
        data = coefficient_label,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        col = "blue"
      )
    }
    p <- p + ggplot2::facet_wrap("label", nrow = 2, scales = "free_y")
  }
  p <- p + ggplot2::theme(
    legend.position = legend_position,
    panel.background = panel_background,
    panel.border = panel_border,
    panel.grid.major = panel_grid_major,
    panel.grid.minor = panel_grid_minor,
    strip.background = strip_background,
  )
  p <- p + ggplot2::xlab(xlab) + ggplot2::ylab(ylab)
  print(p)
  invisible()
}

#' @title Plot the data and the change points for a [fastcpd-class] object
#' @param x A [fastcpd-class] object.
#' @param color_max_count Maximum number of colors to use for the plotting
#' of segments.
#' @param data_point_alpha Alpha of the data points.
#' @param data_point_linewidth Linewidth of the data points.
#' @param data_point_size Size of the data points.
#' @param legend_position Position of the legend.
#' @param panel_background Background of the panel.
#' @param panel_border Border of the panel.
#' @param panel_grid_major Major grid lines of the panel.
#' @param panel_grid_minor Minor grid lines of the panel.
#' @param segment_separator_alpha Alpha of the segment separator lines.
#' @param segment_separator_color Color of the segment separator lines.
#' @param segment_separator_linetype Linetype of the segment separator lines.
#' @param strip_background Background of the strip.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param ... Ignored.
#' @return No return value, called for plotting.
#' @example tests/testthat/examples/plot.R
#'
#' @md
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

#' @title Print the call and the change points for a [fastcpd-class] object
#' @param x A [fastcpd-class] object.
#' @param ... Ignored.
#' @return Return a (temporarily) invisible copy of the [fastcpd-class] object.
#' Called primarily for printing the change points in the model.
#'
#' @md
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

#' @title Show the available methods for a [fastcpd-class] object
#' @param object A [fastcpd-class] object.
#' @return No return value, called for showing a list of available methods
#' for a [fastcpd-class] object.
#'
#' @md
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
    if (!object@cp_only) {
      if (object@family != "custom") {
        cat("\nCost values:\n")
        cat(object@cost_values, "\n")
      }
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

#' @title Show the summary of a [fastcpd-class] object
#' @param object A [fastcpd-class] object.
#' @param ... Ignored.
#' @return Return a (temporarily) invisible copy of the [fastcpd-class] object.
#' Called primarily for printing the summary of the model including the
#' call, the change points, the cost values and the estimated parameters.
#'
#' @md
#' @rdname summary
#' @export
setMethod("summary", "fastcpd", summary.fastcpd)
