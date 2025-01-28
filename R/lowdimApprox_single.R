
#' PCA approximation
#'
#' @param X The data matrix
#' @param ncomp The number of components to use for PCA
#'
#' @return The low dimensional data matrix using PCA
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom mixOmics pca
PCA_approx <- function(X, ncomp = 1) {
  ## checking which taxa has zero variance to avoid errors in pca##
  zero_var_index <- which(X %>% apply(., 2, var) == 0)

  if (length(zero_var_index) > 0) {
    pca_res <- pca(X[, -zero_var_index], ncomp = ncomp, scale = TRUE)

    # loadings for nonZero variance
    loadings_m <- pca_res$loadings$X

    # Including a loading of zero to taxa with zero variance
    for (i in 1:(length(zero_var_index))) {
      if (zero_var_index[i] == 1) {
        loadings_m <- rbind(rep(0, ncomp), loadings_m)
      } else {
        loadings_m <- rbind(
          loadings_m[1:(zero_var_index[i] - 1), ],
          rep(0, ncomp),
          loadings_m[zero_var_index[i]:dim(loadings_m)[1], ]
        )
      }
    }
  } else {
    pca_res <- pca(X, ncomp = ncomp, scale = TRUE)

    # loading matrix
    loadings_m <- pca_res$loadings$X
  }
  return(loadings_m)
}

#' RPCA approximation
#'
#' @param X The data matrix
#' @param ncomp The number of components to use for RPCA
#'
#' @return The low dimensional data matrix Uusing RPCA
#' @export
#'
RPCA_approx <- function(X, ncomp = 1) {
  robpca_res <- rospca::robpca(X, k = ncomp)

  # loading matrix
  loadings_m <- robpca_res$loadings

  return(loadings_m)
}

#' ICA approximation
#'
#' @param X The data matrix
#' @param ncomp The number of components to use for ICA
#'
#' @return The low dimensional data matrix using ICA
#' @export
#'
ICA_approx <- function(X, ncomp = 1) {
  ica_res <- fastICA::fastICA(X, n.comp = ncomp)

  # loading matrix
  loadings_m <- ica_res$K

  return(loadings_m)
}
