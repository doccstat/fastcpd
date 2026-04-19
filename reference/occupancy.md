# Occupancy Detection Data Set

Data set for binary classification of room occupancy from temperature,
humidity, light and CO2 measurements. Ground-truth occupancy was
obtained from time stamped pictures that were taken every minute.

## Usage

``` r
occupancy
```

## Format

A data frame with 9752 rows and 7 variables:

- date:

  Character in the format "YYYY-MM-DD hh:mm:ss" from 2015-02-11 14:48:00
  to 2015-02-18 09:19:00

- Temperature:

  Temperature in Celsius

- Humidity:

  Humidity

- Light:

  Light

- CO2:

  CO2

- HumidityRatio:

  Humidity Ratio

- Occupancy:

  Binary variable with values 0 (unoccupied) and 1

## Source

\<https://github.com/LuisM78/Occupancy-detection-data\>
