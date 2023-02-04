// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
// #include "RcppArmadillo.h"

// [[Rcpp::export]]
// arma::mat data_gen(int n, int d, arma::mat true_coef, arma::vec true_cp_loc, arma::mat Sigma, std::string family, double evar = 0.0) {
//     arma::vec loc = arma::unique(arma::join_cols(arma::join_cols(arma::zeros<arma::vec>(1), true_cp_loc), arma::zeros<arma::vec>(1) + n));
//     arma::mat m1 = arma::eye<arma::mat>(3, 3);
//     arma::mat m2 = arma::eye<arma::mat>(3, 3);
	                     
//     return m1 + 3 * (m1 + m2);
// }

// data_gen <- function(n, d, true.coef, true.cp.loc, Sigma, family, evar = NULL) {
//   loc <- unique(c(0, true.cp.loc, n))
//   if (dim(true.coef)[2] != length(loc) - 1) {
//     stop("true.coef and true.cp.loc do not match")
//   }
//   x <- mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = Sigma)
//   y <- NULL
//   for (i in 1:(length(loc) - 1)) {
//     Xb <- x[(loc[i] + 1):loc[i + 1], , drop = FALSE] %*% true.coef[, i, drop = FALSE]
//     y_new <- if (family == "binomial") {
//       stats::rbinom(length(Xb), 1, 1 / (1 + exp(-Xb)))
//     } else if (family == "poisson") {
//       stats::rpois(length(Xb), exp(Xb))
//     } else if (family == "gaussian") {
//       stats::rnorm(length(Xb), sd = sqrt(evar)) + Xb
//     } else {
//       stop("family not supported")
//     }
//     y <- c(y, y_new)
//   }
//   data <- cbind(x, y)
//   true_cluster <- rep(1:(length(loc) - 1), diff(loc))
//   result <- list(data, true_cluster)
//   return(result)
// }
