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

#' Transcription Profiling of 57 Human Bladder Carcinoma Samples
#'
#' Transcriptome analysis of 57 bladder carcinomas on Affymetrix HG-U95A and
#' HG-U95Av2 microarrays
#'
#' @format A data frame with 2215 rows and 43 variables:
#' \describe{
#'   \item{columns}{Corresponding to the 43 individuals}
#' }
#' @source <https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-TABM-147>
#' @source <https://github.com/cran/ecp/tree/master/data>
#' @example tests/testthat/examples/transcriptome.R
"transcriptome"

#' Occupancy Detection Data Set
#'
#' Data set for binary classification of room occupancy from temperature,
#' humidity, light and CO2 measurements. Ground-truth occupancy was obtained
#' from time stamped pictures that were taken every minute.
#'
#' @format A data frame with 9752 rows and 7 variables:
#' \describe{
#'   \item{date}{Character in the format "YYYY-MM-DD hh:mm:ss" from
#'     2015-02-11 14:48:00 to 2015-02-18 09:19:00}
#'   \item{Temperature}{Temperature in Celsius}
#'   \item{Humidity}{Humidity}
#'   \item{Light}{Light}
#'   \item{CO2}{CO2}
#'   \item{HumidityRatio}{Humidity Ratio}
#'   \item{Occupancy}{Binary variable with values 0 (unoccupied) and 1}
#' }
#' @source <https://github.com/LuisM78/Occupancy-detection-data>
"occupancy"

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
