#' LUPINE function
#'
#' @param data A 3D array of counts or transformed data (e.g. clr) with dimensions samples x taxa x time points
#' @param is.transformed A logical indicating whether the data is transformed or not
#' @param lib_size A matrix of library sizes for each sample and time point
#' @param ncomp The number of components to use for dimensionality reduction
#' @param single A logical indicating whether to use LUPINE for a single time point or longitudinal data
#' @param singleMethod The method to use for dimensionality reduction for single time point
#' @param cutoff The cutoff value for the p-value to determine significance of the correlation
#' @param excluded_taxa A list of taxa to be excluded at each time point
#'
#'
#'
#' @return A list of correlations and p-values for all days
#' @export
#'
LUPINE <- function(data, is.transformed = FALSE, lib_size = NULL, ncomp = 1,
                   single = FALSE, singleMethod = "pca", excluded_taxa = NULL, cutoff = 0.05) {
  if(length(dim(data)) == 2) {
    stop("Data should be a 3D array with dimensions samples x taxa x time points")
  }

  if(is.null(colnames(data))) {
  stop("Data should have taxa names in columns")
  }

  nDays <- dim(data)[3]

  if (single) {
    res <- sapply(1:nDays, function(d) {
      net <- LUPINE_single(data,
        day = d, excluded_taxa, is.transformed, lib_size,
        method = singleMethod, ncomp
      )$pvalue
      net <- apply(net<cutoff, c(1, 2), function(x) {
        ifelse(is.na(x), 0, x)
      })
      return(net)
    }, simplify = FALSE)
  } else {
    res <- sapply(2:nDays, function(d) {
      net <- LUPINE_longitudinal(data,
        day = d, excluded_taxa, is.transformed,
        lib_size, ncomp
      )$pvalue
      net <- apply(net<cutoff, c(1, 2), function(x) {
        ifelse(is.na(x), 0, x)
      })
      return(net)
    }, simplify = FALSE)
  }

  return(res)
}
