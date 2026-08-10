# LUPINE for longitudinal data

LUPINE for longitudinal data

## Usage

``` r
LUPINE_t(
  data,
  day_index,
  num_lags = 999,
  excluded_taxa = NULL,
  is.transformed = FALSE,
  lib_size = NULL,
  ncomp = 1
)
```

## Arguments

- data:

  A 3D array of counts or transformed data with dimensions samples x
  taxa x time points

- day_index:

  The index of the time point for which the correlations are calculated

- num_lags:

  The number of previous time points to be used for the calculation of
  the correlation

- excluded_taxa:

  A list of taxa to be excluded at each time point

- is.transformed:

  A logical indicating whether the data is transformed or not

- lib_size:

  A matrix of library sizes for each sample and time point

- ncomp:

  The number of components to use for dimensionality reduction

## Value

A list of correlations and p-values
