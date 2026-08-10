# LUPINE_bootsrap function

LUPINE_bootsrap function

## Usage

``` r
LUPINE_bootsrap(
  data,
  day_range = NULL,
  is.transformed = FALSE,
  lib_size = NULL,
  ncomp = 1,
  single = FALSE,
  singleMethod = "pca",
  excluded_taxa = NULL,
  cutoff = 0.05,
  nboot = 1000
)
```

## Arguments

- data:

  The data matrix (rows are samples, columns are features/variables, and
  slices are time points)

- day_range:

  day_range The range of days to perform LUPINE bootstrap on

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

  cutoff The cutoff value for the p-value to determine significance of
  the correlation

- nboot:

  The number of bootstrap iterations to perform

## Value

A list of correlations and p-values for all days
