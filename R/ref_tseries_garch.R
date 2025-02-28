# This is a modified copy of tseries/R/garch.R

## Copyright (C) 1997-1999  Adrian Trapletti
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.

garch <-
  function(
      x, order = c(1, 1), series = NULL, control = garch.control(...), ...) {
    if (NCOL(x) > 1) {
      stop("x is not a vector or univariate time series")
    }
    if (!is.vector(order)) stop("order is not a vector")
    switch(control$grad,
      analytical = (agrad <- TRUE),
      numerical = (agrad <- FALSE)
    )
    if (is.null(series)) series <- deparse(substitute(x))
    ists <- stats::is.ts(x)
    x <- stats::as.ts(x)
    xfreq <- stats::frequency(x)
    if (any(is.na(x))) stop("NAs in x")
    if (ists) xtsp <- stats::tsp(x)
    x <- as.matrix(x)
    n <- nrow(x)
    ncoef <- order[1] + order[2] + 1
    small <- 0.05
    coef <- control$start
    if (is.null(coef)) {
      coef <- c(
        stats::var(x) * (1.0 - small * (ncoef - 1)),
        rep.int(small, ncoef - 1)
      )
    }
    if (!is.vector(coef)) stop("coef is not a vector")
    if (ncoef != length(coef)) stop("incorrect length of coef")
    fit <- fit_garch_wrapper(
      as.vector(x, mode = "double"),
      as.integer(n),
      as.vector(coef, mode = "double"),
      as.integer(order[1]),
      as.integer(order[2]),
      as.integer(control$maxiter),
      as.double(control$abstol),
      as.double(control$reltol),
      as.double(control$xtol),
      as.double(control$falsetol),
      as.integer(agrad),
      as.integer(control$trace)
    )
    pred <- pred_garch_wrapper(
      as.vector(x, mode = "double"),
      as.integer(n),
      as.vector(fit$coef, mode = "double"),
      as.integer(order[1]),
      as.integer(order[2]),
      as.integer(FALSE)
    )
    com.hess <- ophess_garch_wrapper(
      as.vector(x, mode = "double"),
      as.integer(n),
      as.vector(fit$coef, mode = "double"),
      as.integer(order[1]),
      as.integer(order[2])
    )
    rank <- do.call(qr, c(list(x = com.hess), control$qr))$rank
    if (rank != ncoef) {
      vc <- matrix(NA, nrow = ncoef, ncol = ncoef)
      warning("singular information")
    } else {
      vc <- solve(com.hess)
    }
    sigt <- sqrt(pred)
    sigt[1:max(order[1], order[2])] <- rep.int(NA, max(order[1], order[2]))
    f <- cbind(sigt, -sigt)
    colnames(f) <- c("sigt", "-sigt")
    e <- as.vector(x) / sigt
    if (ists) {
      attr(e, "tsp") <- attr(f, "tsp") <- xtsp
      attr(e, "class") <- attr(f, "class") <- "ts"
    }
    names(order) <- c("p", "q")
    coef <- fit$coef
    nam.coef <- "a0"
    if (order[2] > 0) {
      nam.coef <- c(nam.coef, paste("a", seq(order[2]), sep = ""))
    }
    if (order[1] > 0) {
      nam.coef <- c(nam.coef, paste("b", seq(order[1]), sep = ""))
    }
    names(coef) <- nam.coef
    colnames(vc) <- rownames(vc) <- nam.coef
    garch <- list(
      order = order,
      coef = coef,
      n.likeli = fit$nlikeli,
      n.used = n,
      residuals = e,
      fitted.values = f,
      series = series,
      frequency = xfreq,
      call = match.call(),
      vcov = vc
    )
    class(garch) <- "garch"
    return(garch)
  }

garch.control <-
  function(
      maxiter = 200, trace = TRUE, start = NULL, grad = c("analytical", "numerical"),
      abstol = max(1e-20, .Machine$double.eps^2),
      reltol = max(1e-10, .Machine$double.eps^(2 / 3)),
      xtol = sqrt(.Machine$double.eps),
      falsetol = 1e2 * .Machine$double.eps, ...) {
    rval <- list(
      maxiter = maxiter, trace = trace, start = start, grad = match.arg(grad),
      abstol = abstol, reltol = reltol, xtol = xtol, falsetol = falsetol
    )
    rval$qr <- list(...)
    rval
  }

# summary.garch <-
# function(object, ...)
# {
#     if(!inherits(object, "garch"))
#         stop("method is only for garch objects")
#     ans <- NULL
#     ans$residuals <- na.remove(object$residuals)
#     tval <- object$coef / sqrt(diag(object$vcov))
#     ans$coef <- cbind(object$coef, sqrt(diag(object$vcov)), tval,
#                       2*(1-pnorm(abs(tval))))
#     dimnames(ans$coef) <-
#         list(names(object$coef),
#              c(" Estimate"," Std. Error"," t value","Pr(>|t|)"))
#     ans$call <- object$call
#     ans$order <- object$order
#     Residuals <- ans$residuals
#     ans$j.b.test <- jarque.bera.test(Residuals)
#     Squared.Residuals <- ans$residuals^2
#     ans$l.b.test <- Box.test(Squared.Residuals, type = "Ljung-Box")
#     class(ans) <- "summary.garch"
#     return(ans)
# }

# logLik.garch <-
# function(object, ...)
# {
#     if(!inherits(object, "garch"))
#         stop("method is only for garch objects")
#     n <- length(na.remove(object$residuals))
#     val <- (-object$n.likeli) - 0.5*n*log(2*pi)
#     attr(val, "df") <- length(object$coef)
#     class(val) <- "logLik"
#     return(val)
# }
