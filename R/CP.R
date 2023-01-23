#' Dynamic Programming with Pruning and Gradient Descent
#'
#' @param data TODO
#' @param beta TODO
#' @param B TODO
#' @param trim TODO
#' @param family TODO
#' @param ... TODO
#'
#' @return TODO
#' @export
CP <- function(data, beta, B = 10, trim = 0.025, family, ...) {
  args_list <- list(...)
  n <- dim(data)[1]
  p <- dim(data)[2] - 1
  Fobj <- c(-beta, 0)
  cp_set <- list(NULL, 0)
  set <- c(0, 1)

  # choose the initial values based on pre-segmentation

  index <- rep(1:B, rep(n / B, B))
  coef.int <- matrix(NA, B, p)
  for (i in 1:B) {
    out <- fastglm::fastglm(
      as.matrix(data[index == i, 1:p]),
      data[index == i, p + 1],
      family
    )
    coef.int[i, ] <- coef(out)
  }
  X1 <- data[1, 1:p]
  if (family == "binomial") {
    cum_coef <- coef <- matrix(coef.int[1, ], p, 1)
    e_eta <- exp(coef %*% X1)
    const <- e_eta / (1 + e_eta)^2
  } else if (family == "poisson") {
    cum_coef <- coef <- DescTools::Winsorize(matrix(coef.int[1, ], p, 1), minval = args_list$L, maxval = args_list$H)
    e_eta <- exp(coef %*% X1)
    const <- e_eta
  }
  cmatrix <- array((X1 %o% X1) * as.numeric(const), c(p, p, 1))

  for (t in 2:n) {
    m <- length(set)
    cval <- rep(NA, m)

    for (i in 1:(m - 1)) {
      coef_c <- coef[, i]
      cum_coef_c <- cum_coef[, i]
      cmatrix_c <- cmatrix[, , i]
      if (family == "binomial") {
        out <- cost_logistic_update(data[t, ], coef_c, cum_coef_c, cmatrix_c)
        coef[, i] <- out[[1]]
        cum_coef[, i] <- out[[2]]
        cmatrix[, , i] <- out[[3]]
        k <- set[i] + 1
        cval[i] <- 0
        if (t - k >= p - 1) {
          cval[i] <- neg_log_lik(data[k:t, ], cum_coef[, i] / (t - k + 1), family = family)
        }
      } else if (family == "poisson") {
        out <- cost_poisson_update(data[t, ], coef_c, cum_coef_c, cmatrix_c, epsilon = args_list$epsilon, G = args_list$G, L = args_list$L, H = args_list$H)
        coef[, i] <- out[[1]]
        cum_coef[, i] <- out[[2]]
        cmatrix[, , i] <- out[[3]]
        k <- set[i] + 1
        cum_coef_win <- DescTools::Winsorize(cum_coef[, i] / (t - k + 1), minval = args_list$L, maxval = args_list$H)
        cval[i] <- 0
        if (t - k >= p - 1) {
          cval[i] <- neg_log_lik(data[k:t, ], cum_coef_win, family = family)
        }
      }
    }

    # the choice of initial values requires further investigation

    cval[m] <- 0
    Xt <- data[t, 1:p]
    if (family == "binomial") {
      cum_coef_add <- coef_add <- coef.int[index[t], ]
      e_eta_t <- exp(coef_add %*% Xt)
      const <- e_eta_t / (1 + e_eta_t)^2
    } else if (family == "poisson") {
      cum_coef_add <- coef_add <- DescTools::Winsorize(coef.int[index[t], ], minval = args_list$L, maxval = args_list$H) ####
      e_eta_t <- exp(coef_add %*% Xt)
      const <- e_eta_t
    }

    cmatrix_add <- (Xt %o% Xt) * as.numeric(const)

    coef <- cbind(coef, coef_add)
    cum_coef <- cbind(cum_coef, cum_coef_add)
    cmatrix <- abind::abind(cmatrix, cmatrix_add, along = 3)

    # Adding a momentum term (TBD)

    obj <- cval + Fobj[set + 1] + beta
    min_val <- min(obj)
    ind <- which(obj == min_val)[1]
    cp_set_add <- c(cp_set[[set[ind] + 1]], set[ind])
    cp_set <- append(cp_set, list(cp_set_add))
    ind2 <- (cval + Fobj[set + 1]) <= min_val
    set <- c(set[ind2], t)
    coef <- coef[, ind2, drop = FALSE]
    cum_coef <- cum_coef[, ind2, drop = FALSE]
    cmatrix <- cmatrix[, , ind2, drop = FALSE]
    Fobj <- c(Fobj, min_val)
  }

  # Remove change-points close to the boundaries

  cp <- cp_set[[n + 1]]
  if (family == "binomial") {
    if (length(cp) > 0) {
      ind3 <- (seq_len(length(cp)))[(cp < trim * n) | (cp > (1 - trim) * n)]
      cp <- cp[-ind3]
    }

    nLL <- 0
    cp_loc <- unique(c(0, cp, n))
    for (i in 1:(length(cp_loc) - 1))
    {
      seg <- (cp_loc[i] + 1):cp_loc[i + 1]
      data_seg <- data[seg, ]
      out <- fastglm::fastglm(
        as.matrix(data_seg[, 1:p]),
        data_seg[, p + 1],
        family
      )
      nLL <- out$deviance / 2 + nLL
    }

    output <- list(cp, nLL)
    names(output) <- c("cp", "nLL")
  } else if (family == "poisson") {
    if (length(cp) > 0) {
      ind3 <- (1:length(cp))[(cp < trim * n) | (cp > (1 - trim) * n)]
      if (length(ind3) > 0) cp <- cp[-ind3]
    }

    cp <- sort(unique(c(0, cp)))
    index <- which((diff(cp) < trim * n) == TRUE)
    if (length(index) > 0) cp <- floor((cp[-(index + 1)] + cp[-index]) / 2)
    cp <- cp[cp > 0]

    # nLL <- 0
    # cp_loc <- unique(c(0,cp,n))
    # for(i in 1:(length(cp_loc)-1))
    # {
    #   seg <- (cp_loc[i]+1):cp_loc[i+1]
    #   data_seg <- data[seg,]
    #   out <- fastglm(as.matrix(data_seg[, 1:p]), data_seg[, p+1], family="Poisson")
    #   nLL <- out$deviance/2 + nLL
    # }

    # output <- list(cp, nLL)
    # names(output) <- c("cp", "nLL")

    output <- list(cp)
    names(output) <- c("cp")
  }

  return(output)
}
