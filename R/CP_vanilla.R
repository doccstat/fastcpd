#' Dynamic Programming with Pruning
#'
#' @param data TODO
#' @param beta TODO
#' @param cost TODO
#'
#' @return TODO
#' @export
CP_vanilla <- function(data, beta, cost = cost_glm) {
  n <- dim(data)[1]
  p <- dim(data)[2] - 1
  Fobj <- c(-beta, 0)
  cp_set <- list(NULL, 0)
  set <- c(0, 1)
  for (t in 2:n) {
    m <- length(set)
    cval <- rep(NA, m)
    for (i in 1:m) {
      k <- set[i] + 1
      cval[i] <- if (t - k >= p - 1) suppressWarnings(cost(data[k:t, ])) else 0
    }
    obj <- cval + Fobj[set + 1] + beta
    min_val <- min(obj)
    ind <- which(obj == min_val)[1]
    cp_set_add <- c(cp_set[[set[ind] + 1]], set[ind])
    cp_set <- append(cp_set, list(cp_set_add))
    ind2 <- (cval + Fobj[set + 1]) <= min_val
    set <- c(set[ind2], t)
    Fobj <- c(Fobj, min_val)
    if (t %% 100 == 0) print(t)
  }
  cp <- cp_set[[n + 1]]
  nLL <- 0
  cp_loc <- unique(c(0, cp, n))
  for (i in 1:(length(cp_loc) - 1)) {
    seg <- (cp_loc[i] + 1):cp_loc[i + 1]
    data_seg <- data[seg, ]
    x <- as.matrix(data_seg[, 1:p])
    out <- fastglm::fastglm(x = x, y = data_seg[, p + 1], family = "binomial")
    nLL <- out$deviance / 2 + nLL
  }

  output <- list(cp, nLL)
  names(output) <- c("cp", "nLL")
  return(output)
}
