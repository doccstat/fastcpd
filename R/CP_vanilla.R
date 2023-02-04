#' Dynamic Programming with Pruning
#'
#' @param data TODO
#' @param beta TODO
#' @param family TODO
#' @param ... TODO
#'
#' @return TODO
#' @export
CP_vanilla <- function(data, beta, family, ...) {
  args_list <- list(...)
  n <- nrow(data)
  p <- ncol(data) - 1

  if (family == "gaussian") {
    B <- args_list$B
    index <- rep(1:B, rep(n / B, B))
    err_sd <- act_num <- rep(NA, B)
    for (i in 1:B) {
      cvfit <- glmnet::cv.glmnet(as.matrix(data[index == i, -1]), data[index == i, 1], family = family)
      coef <- coef(cvfit, s = "lambda.1se")[-1]
      resi <- data[index == i, 1] - as.matrix(data[index == i, -1]) %*% as.numeric(coef)
      err_sd[i] <- sqrt(mean(resi^2))
      act_num[i] <- sum(abs(coef) > 0)
    }
    err_sd_mean <- mean(err_sd) # only works if error sd is unchanged.
    act_num_mean <- mean(act_num)
    beta <- (act_num_mean + 1) * beta # seems to work but there might be better choices
  }

  f_obj <- c(-beta, 0)
  cp_set <- list(NULL, 0)
  set <- c(0, 1)
  for (t in 2:n) {
    m <- length(set)
    cval <- rep(NA, m)
    for (i in 1:m) {
      k <- set[i] + 1
      cval[i] <- 0
      if (family %in% c("binomial", "poisson") && t - k >= p - 1) {
        cval[i] <- suppressWarnings(cost(data[k:t, ], family = family))
      } else if (family == "gaussian" && t - k >= 1) {
        cval[i] <- suppressWarnings(cost(data[k:t, ], family = family, lambda = err_sd_mean * sqrt(2 * log(p) / (t - k + 1))))
      }
    }
    obj <- cval + f_obj[set + 1] + beta
    min_val <- min(obj)
    ind <- which(obj == min_val)[1]
    cp_set_add <- c(cp_set[[set[ind] + 1]], set[ind])
    cp_set <- append(cp_set, list(cp_set_add))
    ind2 <- (cval + f_obj[set + 1]) <= min_val
    set <- c(set[ind2], t)
    f_obj <- c(f_obj, min_val)
    if (t %% 100 == 0) print(t)
  }
  cp <- cp_set[[n + 1]]
  if (family %in% c("binomial", "poisson")) {
    nLL <- 0
    cp_loc <- unique(c(0, cp, n))
    for (i in 1:(length(cp_loc) - 1)) {
      seg <- (cp_loc[i] + 1):cp_loc[i + 1]
      data_seg <- data[seg, ]
      x <- as.matrix(data_seg[, -1])
      out <- fastglm::fastglm(x, data_seg[, 1], family)
      nLL <- out$deviance / 2 + nLL
    }

    output <- list(cp, nLL)
    names(output) <- c("cp", "nLL")
  } else if (family == "gaussian") {
    output <- list(cp)
    names(output) <- c("cp")
  }
  return(output)
}
