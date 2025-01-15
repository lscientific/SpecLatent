## File Name: flatten.R
## File Version: 0.01


#' @title Flattens polytomous matrix
#' @description Flatten a polytomous matrix to a fat binary matrix
#'
#' @param R integer matrix. The polytomous response data matrix. 
#'
#' @return 
#' \describe{
#'   \item{\code{R_flattened}    flattened binary matrix.}
#' }
#' 
#' @export

flatten <- function(R) { 
  N <- nrow(R); J <- ncol(R)
  # categories for each item
  categories <- apply(R, 2, function(x) sort(unique(x)))
  # the number of categories for each item
  Cs <- sapply(categories, length)
  
  R_flattened <- matrix(NA, N, sum(Cs))
  for(i in 1:N) {
    end <- 0
    for(j in 1:J) {
      start <- end + 1
      end <- end + Cs[j]
      
      r <- rep(0, Cs[j])
      r[which(R[i, j] == categories[[j]])] <- 1
      
      R_flattened[i, start:end] <- r
    }
  }
  return(R_flattened)
}

#' @title Re-scale the item parameter matrix 
#' @description Re-scale the item parameter estimation \code{T_hat} for polytomous GoM
#'
#' @param T_mat Numeric matrix. Item parameter matrix.
#' @param Cs Integer vector. The number of categories for each item
#' @return 
#' \describe{
#'   \item{\code{T_mat}    flattened item parameter matrix estimation.}
#' }
#' 
#' @export
rescale_T <- function(T_mat, Cs) {
  J <- nrow(T_mat)
  T_mat[T_mat < 0] <- 0
  
  end <- 0
  for (j in 1:length(Cs)) {
    start <- end + 1
    end <- end + Cs[j]
    mat <- T_mat[start:end, ]
    T_mat[start:end, ] <- apply(mat, 2, function(x) x / sum(x))
  }
  return(T_mat)
}

