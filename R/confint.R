#' Confidence intervals for a fastcpd object
#'
#' @param object A [fastcpd-class] object.
#' @param parm The target parameter. Use \code{"cp"} for change-point
#'   locations or \code{"theta"} for segment parameters.
#' @param level Confidence level.
#' @param method Method used to construct intervals. Change-point intervals
#'   support \code{"bootstrap"} and \code{"profile"}. Parameter intervals
#'   support \code{"wald"}.
#' @param B Number of bootstrap replicates when \code{method = "bootstrap"}.
#' @param bootstrap Bootstrap type. Currently \code{"nonparametric"} resamples
#'   observations within each estimated segment and is available for all
#'   families that can be refitted from \code{object@call}.
#' @param window Optional half-width around each detected change point for
#'   profile intervals. If \code{NULL}, the whole interval between neighboring
#'   detected change points is profiled.
#' @param min_segment_length Minimum number of observations on each side of a
#'   candidate split for profile intervals.
#' @param seed Optional random seed used for bootstrap reproducibility.
#' @param refit_envir Environment used to evaluate bootstrap refits.
#' @param ... Ignored.
#' @return A data frame containing estimates, lower and upper interval bounds,
#'   and method-specific diagnostics.
#'
#' @example tests/testthat/examples/confint.R
#'
#' @md
#' @method confint fastcpd
#' @export
confint.fastcpd <- function(
  object,
  parm = c("cp", "theta"),
  level = 0.95,
  method = NULL,
  B = 999,
  bootstrap = c("nonparametric"),
  window = NULL,
  min_segment_length = 2L,
  seed = NULL,
  refit_envir = parent.frame(),
  ...
) {
  stopifnot("`object` must be a fastcpd object." = methods::is(object, "fastcpd"))
  stopifnot("`level` must be a number in (0, 1)." =
    is.numeric(level) && length(level) == 1 && level > 0 && level < 1)

  parm <- match.arg(parm)
  method_choices <- switch(
    parm,
    cp = c("bootstrap", "profile"),
    theta = "wald"
  )
  if (is.null(method)) {
    method <- method_choices[1]
  }
  method <- match.arg(method, method_choices)

  switch(
    paste(parm, method, sep = ":"),
    "cp:bootstrap" = fastcpd_confint_cp_bootstrap(
      object = object,
      level = level,
      B = B,
      bootstrap = bootstrap,
      seed = seed,
      refit_envir = refit_envir
    ),
    "cp:profile" = fastcpd_confint_cp_profile(
      object = object,
      level = level,
      window = window,
      min_segment_length = min_segment_length
    ),
    "theta:wald" = fastcpd_confint_theta_wald(
      object = object,
      level = level
    )
  )
}

fastcpd_confint_cp_bootstrap <- function(
  object,
  level,
  B,
  bootstrap,
  seed,
  refit_envir
) {
  bootstrap <- match.arg(bootstrap, "nonparametric")
  B <- as.integer(B)
  stopifnot("`B` must be a positive integer." = length(B) == 1 && B > 0)

  cp_set <- sort(as.integer(object@cp_set))
  if (!length(cp_set)) {
    return(fastcpd_empty_confint("cp", level, "bootstrap"))
  }

  if (!is.null(seed)) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
    }
    on.exit({
      if (has_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  matched <- matrix(NA_real_, nrow = B, ncol = length(cp_set))
  failed_refits <- 0L
  for (b in seq_len(B)) {
    boot_data <- fastcpd_segment_bootstrap_data(object, bootstrap)
    boot_result <- tryCatch(
      fastcpd_refit_with_data(object, boot_data, cp_only = TRUE, refit_envir),
      error = function(e) e
    )
    if (!inherits(boot_result, "error")) {
      matched[b, ] <- fastcpd_match_cp_set(
        cp_set,
        sort(as.integer(boot_result@cp_set)),
        nrow(object@data)
      )
    } else {
      failed_refits <- failed_refits + 1L
    }
  }
  if (failed_refits > 0L) {
    warning(
      failed_refits, " bootstrap refit(s) failed; their intervals use `NA`.",
      call. = FALSE
    )
  }

  alpha <- 1 - level
  lower <- upper <- detection_rate <- numeric(length(cp_set))
  for (i in seq_along(cp_set)) {
    estimates <- matched[, i]
    detection_rate[i] <- mean(!is.na(estimates))
    if (all(is.na(estimates))) {
      lower[i] <- NA_real_
      upper[i] <- NA_real_
    } else {
      interval <- stats::quantile(
        estimates,
        probs = c(alpha / 2, 1 - alpha / 2),
        na.rm = TRUE,
        names = FALSE,
        type = 1
      )
      lower[i] <- interval[1]
      upper[i] <- interval[2]
    }
  }

  data.frame(
    parm = "cp",
    index = seq_along(cp_set),
    estimate = cp_set,
    lower = lower,
    upper = upper,
    detection_rate = detection_rate,
    level = level,
    method = "bootstrap",
    bootstrap = bootstrap,
    row.names = NULL
  )
}

fastcpd_confint_cp_profile <- function(
  object,
  level,
  window,
  min_segment_length
) {
  cp_set <- sort(as.integer(object@cp_set))
  if (!length(cp_set)) {
    return(fastcpd_empty_confint("cp", level, "profile"))
  }

  min_segment_length <- as.integer(min_segment_length)
  stopifnot("`min_segment_length` must be a positive integer." =
    length(min_segment_length) == 1 && min_segment_length > 0)

  cost_function <- fastcpd_profile_cost_function(object)
  n <- nrow(object@data)
  bounds <- c(0L, cp_set, n)
  cutoff <- stats::qchisq(level, df = 1) / 2
  out <- vector("list", length(cp_set))

  for (i in seq_along(cp_set)) {
    left <- bounds[i] + 1L
    right <- bounds[i + 2L]
    tau_min <- left + min_segment_length - 1L
    tau_max <- right - min_segment_length
    if (!is.null(window)) {
      tau_min <- max(tau_min, cp_set[i] - as.integer(window))
      tau_max <- min(tau_max, cp_set[i] + as.integer(window))
    }
    candidates <- seq.int(tau_min, tau_max)
    if (!length(candidates)) {
      out[[i]] <- data.frame(
        parm = "cp", index = i, estimate = cp_set[i],
        lower = NA_real_, upper = NA_real_, profile_min = NA_real_,
        cutoff = cutoff, level = level, method = "profile"
      )
      next
    }

    costs <- vapply(
      candidates,
      function(tau) {
        cost_function(left, tau) + cost_function(tau + 1L, right)
      },
      numeric(1)
    )
    finite <- is.finite(costs)
    if (!any(finite)) {
      lower <- upper <- profile_min <- NA_real_
    } else {
      profile_min <- min(costs[finite])
      support <- candidates[finite][costs[finite] - profile_min <= cutoff]
      lower <- min(support)
      upper <- max(support)
    }
    out[[i]] <- data.frame(
      parm = "cp",
      index = i,
      estimate = cp_set[i],
      lower = lower,
      upper = upper,
      profile_min = profile_min,
      cutoff = cutoff,
      level = level,
      method = "profile"
    )
  }

  do.call(rbind, out)
}

fastcpd_confint_theta_wald <- function(object, level) {
  if (isTRUE(object@cp_only) || !nrow(object@thetas) || !ncol(object@thetas)) {
    stop("Wald intervals require a fastcpd object fitted with `cp_only = FALSE`.")
  }

  se_function <- fastcpd_theta_se_function(object)
  theta <- as.matrix(object@thetas)
  bounds <- fastcpd_segment_bounds(object)
  z_value <- stats::qnorm(1 - (1 - level) / 2)
  out <- vector("list", ncol(theta))

  for (segment in seq_len(ncol(theta))) {
    se <- se_function(bounds$start[segment], bounds$end[segment])
    estimate <- theta[, segment]
    out[[segment]] <- data.frame(
      parm = "theta",
      segment = segment,
      parameter = seq_along(estimate),
      estimate = estimate,
      lower = estimate - z_value * se,
      upper = estimate + z_value * se,
      se = se,
      level = level,
      method = "wald"
    )
  }

  do.call(rbind, out)
}

fastcpd_empty_confint <- function(parm, level, method) {
  data.frame(
    parm = parm,
    index = integer(0),
    estimate = numeric(0),
    lower = numeric(0),
    upper = numeric(0),
    level = level,
    method = method
  )
}

fastcpd_segment_bootstrap_data <- function(object, bootstrap) {
  if (bootstrap != "nonparametric") {
    stop("Only `bootstrap = \"nonparametric\"` is currently implemented.")
  }
  data <- object@data
  n <- nrow(data)
  bounds <- c(0L, sort(as.integer(object@cp_set)), n)
  boot_data <- data
  for (i in seq_len(length(bounds) - 1L)) {
    rows <- seq.int(bounds[i] + 1L, bounds[i + 1L])
    if (length(rows)) {
      boot_data[rows, ] <- data[
        sample(rows, length(rows), replace = TRUE),
        ,
        drop = FALSE
      ]
    }
  }
  boot_data
}

fastcpd_refit_with_data <- function(object, data, cp_only, refit_envir) {
  refit_call <- object@call
  if (!is.call(refit_call)) {
    stop("`object@call` is not a refittable call.")
  }

  call_names <- names(refit_call)
  if ("data" %in% call_names) {
    refit_call$data <- data
  } else {
    call_name <- as.character(refit_call[[1L]])[1L]
    data_position <- if (identical(call_name, "fastcpd")) 3L else 2L
    if (length(refit_call) < data_position) {
      stop("Could not locate the data argument in `object@call`.")
    }
    refit_call[[data_position]] <- data
  }
  refit_call$cp_only <- cp_only
  refit_call$show.progress <- FALSE

  eval(refit_call, envir = refit_envir)
}

fastcpd_match_cp_set <- function(reference_cp, bootstrap_cp, n) {
  if (!length(bootstrap_cp)) {
    return(rep(NA_real_, length(reference_cp)))
  }
  matched <- rep(NA_real_, length(reference_cp))
  for (i in seq_along(reference_cp)) {
    left <- if (i == 1L) {
      0
    } else {
      floor((reference_cp[i - 1L] + reference_cp[i]) / 2)
    }
    right <- if (i == length(reference_cp)) {
      n
    } else {
      ceiling((reference_cp[i] + reference_cp[i + 1L]) / 2)
    }
    candidates <- bootstrap_cp[bootstrap_cp > left & bootstrap_cp <= right]
    if (length(candidates)) {
      matched[i] <- candidates[which.min(abs(candidates - reference_cp[i]))]
    }
  }
  matched
}

fastcpd_segment_bounds <- function(object) {
  bounds <- c(0L, sort(as.integer(object@cp_set)), nrow(object@data))
  data.frame(
    segment = seq_len(length(bounds) - 1L),
    start = bounds[-length(bounds)] + 1L,
    end = bounds[-1L]
  )
}

fastcpd_profile_cost_function <- function(object) {
  family <- object@family
  data <- as.matrix(object@data)
  storage.mode(data) <- "double"

  switch(
    family,
    mean = function(start, end) {
      segment <- data[start:end, , drop = FALSE]
      centered <- sweep(segment, 2, colMeans(segment), check.margin = FALSE)
      sum(centered^2) / 2
    },
    variance = {
      centered_data <- sweep(data, 2, colMeans(data), check.margin = FALSE)
      function(start, end) {
        segment <- centered_data[start:end, , drop = FALSE]
        n <- nrow(segment)
        covariance <- crossprod(segment) / n
        n * fastcpd_logdet(covariance) / 2
      }
    },
    meanvariance = function(start, end) {
      segment <- data[start:end, , drop = FALSE]
      n <- nrow(segment)
      centered <- sweep(segment, 2, colMeans(segment), check.margin = FALSE)
      covariance <- crossprod(centered) / n
      n * fastcpd_logdet(covariance) / 2
    },
    exponential = function(start, end) {
      segment <- data[start:end, 1, drop = TRUE]
      if (any(segment <= 0)) return(Inf)
      n <- length(segment)
      n * (log(mean(segment)) + 1)
    },
    lm = fastcpd_profile_cost_lm(data),
    binomial = fastcpd_profile_cost_glm(data, stats::binomial()),
    poisson = fastcpd_profile_cost_glm(data, stats::poisson()),
    quantile = fastcpd_profile_cost_quantile(data, object@order[1]),
    arima = fastcpd_profile_cost_arima(data, object@order),
    stop(sprintf(
      "Profile intervals are not implemented for family `%s`.",
      family
    ))
  )
}

fastcpd_profile_cost_lm <- function(data) {
  function(start, end) {
    segment <- data[start:end, , drop = FALSE]
    y <- segment[, 1]
    x <- segment[, -1, drop = FALSE]
    if (nrow(x) <= ncol(x)) return(Inf)
    fit <- tryCatch(stats::lm.fit(x, y), error = function(e) NULL)
    if (is.null(fit)) return(Inf)
    sum(fit$residuals^2) / 2
  }
}

fastcpd_profile_cost_glm <- function(data, family) {
  function(start, end) {
    segment <- data[start:end, , drop = FALSE]
    y <- segment[, 1]
    x <- segment[, -1, drop = FALSE]
    if (nrow(x) <= ncol(x)) return(Inf)
    fit <- tryCatch(
      suppressWarnings(stats::glm.fit(x = x, y = y, family = family)),
      error = function(e) NULL
    )
    if (is.null(fit) || !is.finite(fit$deviance)) return(Inf)
    fit$deviance / 2
  }
}

fastcpd_profile_cost_quantile <- function(data, tau) {
  function(start, end) {
    segment <- data[start:end, , drop = FALSE]
    y <- segment[, 1]
    x <- segment[, -1, drop = FALSE]
    if (nrow(x) <= ncol(x)) return(Inf)
    fastcpd_quantile_cost(x, y, tau)
  }
}

fastcpd_profile_cost_arima <- function(data, order) {
  function(start, end) {
    segment <- data[start:end, 1, drop = TRUE]
    tryCatch(
      -stats::arima(segment, order = order, method = "ML")$loglik,
      error = function(e) Inf
    )
  }
}

fastcpd_quantile_cost <- function(x, y, tau) {
  beta <- tryCatch(
    solve(crossprod(x) + diag(1e-6, ncol(x)), crossprod(x, y)),
    error = function(e) NULL
  )
  if (is.null(beta)) return(Inf)
  for (i in seq_len(100L)) {
    residual <- c(y - x %*% beta)
    weights <- ifelse(residual > 0, tau, 1 - tau) /
      pmax(abs(residual), 1e-6)
    wx <- x * weights
    beta_new <- tryCatch(
      solve(crossprod(x, wx) + diag(1e-6, ncol(x)), crossprod(wx, y)),
      error = function(e) NULL
    )
    if (is.null(beta_new)) return(Inf)
    delta <- sqrt(sum((beta_new - beta)^2))
    beta <- beta_new
    if (delta < 1e-8 * (1 + sqrt(sum(beta^2)))) break
  }
  residual <- c(y - x %*% beta)
  sum(ifelse(residual >= 0, tau * residual, (tau - 1) * residual))
}

fastcpd_logdet <- function(matrix) {
  if (!all(is.finite(matrix))) {
    return(Inf)
  }
  eigenvalues <- tryCatch(
    eigen(matrix, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )
  if (anyNA(eigenvalues)) return(Inf)
  sum(log(pmax(eigenvalues, .Machine$double.eps)))
}

fastcpd_theta_se_function <- function(object) {
  family <- object@family
  data <- as.matrix(object@data)
  storage.mode(data) <- "double"

  switch(
    family,
    mean = function(start, end) {
      segment <- data[start:end, , drop = FALSE]
      n <- nrow(segment)
      if (n <= 1L) return(rep(NA_real_, ncol(segment)))
      variance <- stats::var(segment)
      if (length(variance) == 1L) {
        sqrt(as.numeric(variance) / n)
      } else {
        sqrt(diag(variance) / n)
      }
    },
    exponential = function(start, end) {
      segment <- data[start:end, 1, drop = TRUE]
      rate <- 1 / mean(segment)
      rate / sqrt(length(segment))
    },
    lm = fastcpd_theta_se_lm(data),
    binomial = fastcpd_theta_se_glm(data, stats::binomial()),
    poisson = fastcpd_theta_se_glm(data, stats::poisson()),
    stop(sprintf(
      "Wald intervals are not implemented for family `%s`.",
      family
    ))
  )
}

fastcpd_theta_se_lm <- function(data) {
  function(start, end) {
    segment <- data[start:end, , drop = FALSE]
    y <- segment[, 1]
    x <- segment[, -1, drop = FALSE]
    fit <- tryCatch(stats::lm.fit(x, y), error = function(e) NULL)
    if (is.null(fit)) return(rep(NA_real_, ncol(x)))
    df <- max(nrow(x) - ncol(x), 1L)
    sigma2 <- sum(fit$residuals^2) / df
    xtx_inv <- tryCatch(solve(crossprod(x)), error = function(e) NULL)
    if (is.null(xtx_inv)) return(rep(NA_real_, ncol(x)))
    sqrt(diag(xtx_inv) * sigma2)
  }
}

fastcpd_theta_se_glm <- function(data, family) {
  function(start, end) {
    segment <- data[start:end, , drop = FALSE]
    y <- segment[, 1]
    x <- segment[, -1, drop = FALSE]
    fit <- tryCatch(
      suppressWarnings(stats::glm.fit(x = x, y = y, family = family)),
      error = function(e) NULL
    )
    if (is.null(fit)) return(rep(NA_real_, ncol(x)))
    xtwx_inv <- tryCatch(
      solve(crossprod(x, x * fit$weights)),
      error = function(e) NULL
    )
    if (is.null(xtwx_inv)) return(rep(NA_real_, ncol(x)))
    sqrt(diag(xtwx_inv))
  }
}
