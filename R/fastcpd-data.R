#' Bitcoin Market Price (USD)
#'
#' The average USD market price across major bitcoin exchanges.
#'
#' @format  A data frame with 1354 rows and 2 variables:
#' \describe{
#'   \item{date}{POSIXct,POSIXt (TZ: "UTC") from 2019-01-02 to 2023-10-28}
#'   \item{price}{The average USD market price across major bitcoin exchanges}
#' }
#' @source <https://www.blockchain.com/explorer/charts/market-price>
#' @example tests/testthat/examples/bitcoin.txt
"bitcoin"

#' Well-log Dataset from Numerical Bayesian Methods Applied to Signal Processing
#'
#' This is the well-known well-log dataset used in many changepoint papers
#' obtained from Alan Turing Institute GitHub repository and licensed under
#' the MIT license. Outliers with value less or equal to 1e5 are removed.
#'
#' @format A Time-Series of length 3989.
#' @source <https://github.com/alan-turing-institute/TCPD>
#' @example tests/testthat/examples/well_log.txt
"well_log"
