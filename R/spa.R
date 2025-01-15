## File Name: spa.R
## File Version: 0.01


#' @title Sequential projection algorithm
#'
#' @description A sequential projection algorithm (SPA) to find the pure subjects
#'
#' @param mat the (pruned) left singular matrix to conduct SPA on. 

#' @return 
#' \describe{
#'   \item{\code{S_hat}    a vector of the pure subject indices.}
#' }
#' 
#' @references Gillis, N. and Vavasis, S. A. (2013). Fast and robust recursive algorithms for separable nonnegative matrix factorization. IEEE Transactions on Pattern Analysis and Machine Intelligence, 36(4):698–714.


spa <- function(mat) {
  K <- ncol(mat)
  S_hat <- rep(NA, K)
  
  Y <- mat  # make a copy of the matrix
  for (k in 1:K) {
    row_norms <- apply(Y, 1, function(x) sqrt(sum(x^2)))  # calculate the row norms
    idx <- which.max(row_norms)  # identify the largest norm
    S_hat[k] <- idx
    u <- Y[idx, ] / row_norms[idx]
    Y <- Y %*% (diag(K) - matrix(u, nc=1) %*% matrix(u, nr=1))  # projection onto the orthogonal subspace
  }
  
  return(S_hat)
}
