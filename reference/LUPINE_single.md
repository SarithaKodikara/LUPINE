# LUPINE for single time point

LUPINE for single time point

## Usage

``` r
LUPINE_single(
  data,
  day,
  excluded_taxa = NULL,
  is.transformed = FALSE,
  lib_size = NULL,
  method = "pca",
  ncomp = 1
)
```

## Arguments

- data:

  A 3D array of counts or transformed data with dimensions samples x
  taxa x time points

- day:

  The time point for which the correlations are calculated

- excluded_taxa:

  A list of taxa to be excluded at each time point

- is.transformed:

  A logical indicating whether the data is transformed or not

- lib_size:

  A matrix of library sizes for each sample and time point

- method:

  The method to use for dimensionality reduction. Options are "pca",
  "ica", "rpca"

- ncomp:

  The number of components to use for dimensionality reduction

## Value

A list of correlations and p-values
