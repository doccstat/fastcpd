% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastcpd.R
\name{fastcpd_meanvariance}
\alias{fastcpd_meanvariance}
\alias{fastcpd.meanvariance}
\alias{fastcpd_mv}
\alias{fastcpd.mv}
\title{Find change points efficiently in mean variance change models}
\usage{
fastcpd_meanvariance(data, ...)

fastcpd.meanvariance(data, ...)

fastcpd_mv(data, ...)

fastcpd.mv(data, ...)
}
\arguments{
\item{data}{A matrix, a data frame or a vector.}

\item{...}{Other arguments passed to \code{\link[=fastcpd]{fastcpd()}}, for example,
\code{segment_count}.}
}
\value{
A \linkS4class{fastcpd} object.
}
\description{
\code{\link[=fastcpd_meanvariance]{fastcpd_meanvariance()}}, \code{\link[=fastcpd.meanvariance]{fastcpd.meanvariance()}},
\code{\link[=fastcpd_mv]{fastcpd_mv()}}, \code{\link[=fastcpd.mv]{fastcpd.mv()}} are wrapper
functions of \code{\link[=fastcpd]{fastcpd()}} to find the meanvariance change. The
function is similar to \code{\link[=fastcpd]{fastcpd()}} except that the data is by
default a matrix or data frame or a vector with each row / element as an
observation and thus a formula is not required here.
}
\examples{
set.seed(1)
data <- c(
  rnorm(3000, 0, 1),
  rnorm(1000, 10, 1),
  rnorm(3000, 10, 20),
  rnorm(1000, 0, 1)
)
system.time(result <- fastcpd.mv(data))
summary(result)
plot(result)
set.seed(1)
p <- 3
data <- if (requireNamespace("mvtnorm", quietly = TRUE)) {
  rbind(
    mvtnorm::rmvnorm(2e+5, mean = rep(0, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(1e+5, mean = rep(50, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(2e+5, mean = rep(0, p), sigma = diag(100, p)),
    mvtnorm::rmvnorm(2e+5, mean = rep(0, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(1e+5, mean = rep(50, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(2e+5, mean = rep(50, p), sigma = diag(100, p))
  )
} else {
  rbind(
    matrix(rnorm(p * 2e+5, mean = 0, sd = 1), ncol = p),
    matrix(rnorm(p * 1e+5, mean = 50, sd = 1), ncol = p),
    matrix(rnorm(p * 2e+5, mean = 0, sd = 10), ncol = p),
    matrix(rnorm(p * 2e+5, mean = 0, sd = 1), ncol = p),
    matrix(rnorm(p * 1e+5, mean = 50, sd = 1), ncol = p),
    matrix(rnorm(p * 2e+5, mean = 50, sd = 10), ncol = p)
  )
}
system.time(result <- fastcpd.mv(data))
summary(result)
result@thetas[seq_len(p), ]
lapply(result@thetas[seq_len(p^2) + p, ], function(thetas) matrix(thetas, p))
}
\seealso{
\code{\link[=fastcpd]{fastcpd()}}
}
