## File Name: perm.R
## File Version: 0.01


#' @title Permutes a vector
#' @description This function performs permutation on a vector
#'
#' @param x Numeric vector. Vector of labels with integer values 1, ..., K
#' @param p Numeric vector. An integer permutation vector.
#' @return Permuted vector \code{x_perm}

perm <- function(x, p) {
  K <- length(p)
  x_perm <- rep(NA, length(x))
  for (i in 1:length(x)) {
    for (k in 1:K) {
      x_perm[x == k] <- p[k]
    }
  }
  return(x_perm)
}
