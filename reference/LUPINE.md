# LUPINE function

LUPINE function

## Usage

``` r
LUPINE(
  data,
  is.transformed = FALSE,
  lib_size = NULL,
  ncomp = 1,
  single = FALSE,
  singleMethod = "pca",
  excluded_taxa = NULL,
  cutoff = 0.05
)
```

## Arguments

- data:

  A 3D array of counts or transformed data (e.g. clr) with dimensions
  samples x taxa x time points

- is.transformed:

  A logical indicating whether the data is transformed or not

- lib_size:

  A matrix of library sizes for each sample and time point

- ncomp:

  The number of components to use for dimensionality reduction

- single:

  A logical indicating whether to use LUPINE for a single time point or
  longitudinal data

- singleMethod:

  The method to use for dimensionality reduction for single time point

- excluded_taxa:

  A list of taxa to be excluded at each time point

- cutoff:

  The cutoff value for the p-value to determine significance of the
  correlation

## Value

A list of correlations and p-values for all days
