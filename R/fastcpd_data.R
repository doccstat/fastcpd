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
#' @example tests/testthat/examples/data-bitcoin.txt
"bitcoin"

#' Transcription Profiling of 57 Human Bladder Carcinoma Samples
#'
#' Transcriptome analysis of 57 bladder carcinomas on Affymetrix HG-U95A and
#' HG-U95Av2 microarrays
#'
#' @format A data frame with 2215 rows and 43 variables:
#' \describe{
#'   \item{3}{Individual 3}
#'   \item{4}{Individual 4}
#'   \item{5}{Individual 5}
#'   \item{6}{Individual 6}
#'   \item{7}{Individual 7}
#'   \item{8}{Individual 8}
#'   \item{9}{Individual 9}
#'   \item{10}{Individual 10}
#'   \item{14}{Individual 14}
#'   \item{15}{Individual 15}
#'   \item{16}{Individual 16}
#'   \item{17}{Individual 17}
#'   \item{18}{Individual 18}
#'   \item{19}{Individual 19}
#'   \item{21}{Individual 21}
#'   \item{22}{Individual 22}
#'   \item{24}{Individual 24}
#'   \item{26}{Individual 26}
#'   \item{28}{Individual 28}
#'   \item{30}{Individual 30}
#'   \item{31}{Individual 31}
#'   \item{33}{Individual 33}
#'   \item{34}{Individual 34}
#'   \item{35}{Individual 35}
#'   \item{36}{Individual 36}
#'   \item{37}{Individual 37}
#'   \item{38}{Individual 38}
#'   \item{39}{Individual 39}
#'   \item{40}{Individual 40}
#'   \item{41}{Individual 41}
#'   \item{42}{Individual 42}
#'   \item{43}{Individual 43}
#'   \item{44}{Individual 44}
#'   \item{45}{Individual 45}
#'   \item{46}{Individual 46}
#'   \item{47}{Individual 47}
#'   \item{48}{Individual 48}
#'   \item{49}{Individual 49}
#'   \item{50}{Individual 50}
#'   \item{51}{Individual 51}
#'   \item{53}{Individual 53}
#'   \item{54}{Individual 54}
#'   \item{57}{Individual 57}
#' }
#' @source <https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-TABM-147>
#' @source <https://github.com/cran/ecp/tree/master/data>
#' @example tests/testthat/examples/data-transcriptome.txt
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

#' UK Seatbelts Data
#'
#' Road Casualties in Great Britain 1969–84.
#'
#' @format \code{uk_seatbelts} is a multiple time series, with columns
#' \describe{
#'  \item{DriversKilled}{car drivers killed.}
#'  \item{front}{front-seat passengers killed or seriously injured.}
#'  \item{rear}{rear-seat passengers killed or seriously injured.}
#'  \item{kms}{distance driven.}
#'  \item{PetrolPrice}{petrol price.}
#'  \item{VanKilled}{number of van (‘light goods vehicle’) drivers.}
#'  \item{law}{0/1: was the law in effect that month?}
#' }
#' @source R package \pkg{datasets}
#' @example tests/testthat/examples/data-uk_seatbelts.R
"uk_seatbelts"

#' Well-log Dataset from Numerical Bayesian Methods Applied to Signal Processing
#'
#' This is the well-known well-log dataset used in many changepoint papers
#' obtained from Alan Turing Institute GitHub repository and licensed under
#' the MIT license. Outliers with value less or equal to 1e5 are removed.
#'
#' @format A Time-Series of length 3989.
#' @source <https://github.com/alan-turing-institute/TCPD>
#' @example tests/testthat/examples/data-well_log_1.R
#' @example tests/testthat/examples/data-well_log_2.R
#' @example tests/testthat/examples/data-well_log-quantile.txt
"well_log"
