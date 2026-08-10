# Example Dataset: HFHS data

The example data is a list containing 8 elements:

1.  OTU data for HFHS diet: an array with dimensions n=23, p=212, t=4.

2.  OTU data for Normal diet: an array with dimensions n=23, p=212, t=4.

3.  Library size for HFHS diet: a matrix with dimensions n=23, t=4.

4.  Library size for Normal diet: a matrix with dimensions n=23, t=4.

5.  Filtered sample information: a data frame with dimensions.

6.  Filtered taxonomy data: a data frame with dimension.

7.  Low abundant taxa for HFHS diet: a list of 4 elements.

8.  Low abundant taxa for Normal diet: a list of 4 elements.

## Usage

``` r
HFHSdata
```

## Format

A list containing 8 elements, including arrays, matrices, data frames,
and lists, as described above.

## References

Kodikara, S., & Le Cao, K. A. (2024). Microbial network inference for
longitudinal microbiome studies with LUPINE. bioRxiv, 2024-05.

## Examples

``` r
data(HFHSdata)
```
