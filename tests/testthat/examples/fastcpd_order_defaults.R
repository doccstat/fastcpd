# Omitted orders select the same model through either public entry point.
default_order_cases <- list(
  ar = 1, arma = c(1, 1), arima = c(1, 1, 0), garch = c(1, 1),
  var = 1, quantile = 0.5, kcp = c(100, 0)
)
default_order_series <- sin(seq_len(24)) + 0.2 * cos(seq_len(24) * 2)
default_order_results <- lapply(names(default_order_cases), function(family) {
  data <- matrix(default_order_series, ncol = 1)
  if (family == "var") data <- cbind(data, cos(seq_len(24) * 0.7))
  if (family == "quantile") data <- cbind(data, 1)
  wrapper <- get(paste0("detect_", family), asNamespace("fastcpd"))
  arguments <- list(data = data, beta = 1e6)
  # Reuse the same random features for KCP's generic and named calls.
  set.seed(31)
  named <- do.call(wrapper, arguments)
  generic_arguments <- c(
    list(formula = ~ . - 1, family = family),
    list(data = as.data.frame(data), beta = 1e6)
  )
  set.seed(31)
  generic <- do.call(detect, generic_arguments)
  set.seed(31)
  explicit <- do.call(detect, c(
    generic_arguments, list(order = default_order_cases[[family]])
  ))
  invalid_order <- rep(0, length(default_order_cases[[family]]))
  invalid_generic <- if (family != "kcp") tryCatch(
    do.call(detect, c(generic_arguments, list(order = invalid_order))),
    error = identity
  ) else NULL
  invalid_named <- if (family != "kcp") tryCatch(
    do.call(wrapper, c(arguments, list(order = invalid_order))),
    error = identity
  ) else NULL
  list(
    named = named, generic = generic, explicit = explicit,
    invalid_generic = invalid_generic, invalid_named = invalid_named
  )
})
names(default_order_results) <- names(default_order_cases)
default_order_mean <- detect_mean(default_order_series, beta = 1e6)
