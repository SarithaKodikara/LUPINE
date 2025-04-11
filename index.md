# LongitUdinal modelling with Partial least squares regression for NEtwork inference (LUPINE)

[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://opensource.org/licenses/GPL-3.0)
[![Project Status: Active](https://img.shields.io/badge/project%20status-active-brightgreen.svg)](https://example.com) 


## Installation

Install the GitHub version with:

```r
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("https://github.com/SarithaKodikara/LUPINE")
```

## What is LUPINE?

LUPINE is an R package designed for inferring associations between two microbial species using partial correlations for longitudinal data. 
It is particularly useful for analysing data from longitudinal studies, where repeated measurements are taken over time. 

## When to use LUPINE?

LUPINE supports two modes: a single time point analysis (`single = TRUE` in `LUPINE()`) and a longitudinal mode that incorporates past data to infer associations. 

The longitudinal mode requires
**matching samples across all time points**. 
Therefore, it cannot be used with missing samples unless you either **impute the missing data** (see [Case Study 3](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-025-02041-w))
 or **remove samples with incomplete data**.
 
As LUPINE is a correlation-based method, **sample size is critical** for obtaining reliable results. 
As a rule of thumb, ensure a **minimum of 5 samples** is available for analysis.

## How to use LUPINE?
You need to load the package and then run the LUPINE function with your data. 
Make sure you read following sections to understand the arguments in the function.


- ### data

LUPINE is designed to work with **3D arrays**, where the first dimension (rows) represents samples, the second dimension (columns) represents microbial species, and the third dimension represents time points.
This array can contain **count data** or **transformed/normalized data** (e.g., log-transformed, CLR-transformed, etc.).
By default, LUPINE assumes count data. If your data is *already transformed*, set the argument `is.transformed = TRUE`.

<code style="color : Blue">  Note </code> 
Your data should not contain taxa with zero variance at any time point.
If any taxa have zero variance, make sure to exclude them in modelling by including them in the `excluded_taxa` argument.

- ### day

A **numeric vector** indicating the time points (e.g., days or weeks) corresponding to each slice along the third dimension of the data array.
The values must be in **increasing order** and should match the order of time points in the array.

- ### excluded_taxa (optional)

A **list of vectors** specifying microbial species to exclude at each time point.
Each vector corresponds to a time point, and the taxa names must **match the column names** of the data array.

This is useful for excluding taxa with low counts or many zeros, which may lead to unstable or biased correlation estimates.

- ### lib_size (optional)

If your input data is counts and you have access to library sizes, 
you can provide them via the lib_size argument as a **matrix** where rows represent samples and columns represent time points.

Ensure the matrix contains **no zero values**, as zero library sizes are not valid for the analysis.
