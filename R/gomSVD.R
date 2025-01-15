## File Name: gomSVD.R
## File Version: 0.01


#' @title gGoM estimation based on the left singular matrix
#'
#' @description Estimation algorithm for gGoM model with the left singular matrix.
#'
#' @param U the pruned left singular matrix from data SVD. 
#' @param V the right singular matrix from data SVD. 
#' @param d the vector containing the singular values. 
#' @param prune logical; if true, the pruning step is performed.
#' @param r the number of neighbors to consider in pruning. Used only when \code{prune} is \code{TRUE}. Default value is 10.
#' @param q the cutoff for the upper quantile of row norms. Used only when \code{prune} is \code{TRUE}. Higher \code{q} leads to more points being pruned. Default value is 0.4.
#' @param e the cutoff for the upper quantile of average distance. Used only when \code{prune} is \code{TRUE}. Higher \code{e} leads to more points being pruned. Default value is 0.2.
#' @return The function returns a list with the following components:
#' \itemize{
#' \item{\code{P_hat} the estimated membership scores.}
#' \item{\code{T_hat} the estimated item response parameters (not truncated).}
#' \item{\code{R_hat} the estimated response expectation (not truncated).}
#' \item{\code{S_hat} the estimated indices of pure subjects.}
#' \item{\code{t} computation time.}
#' }
#' @export

gomSVD <- function(U, V, d, prune=T, r=10, q=0.4, e=0.2) {
  N <- nrow(U)
  if (ncol(U) < 2) stop("Error: the column number of U needs to be at least 2!")
  
  if (prune) {
    indices <- pruning(U, r, q, e)  # obtain the subjects to be pruned
    X <- U[-indices,]
  } else {X <- U}
  
  S_hat_X <- spa(X)  # SPA finds the pure subjects from the pruned matrix
  vertices <- X[S_hat_X, ]  # the simplex vertices
  S_hat <- (1:N)[-indices][S_hat_X]
  
  # estimation for Pi
  P1 <- U %*% solve(vertices)  
  P2 <- t(apply(P1, 1, function(x) ifelse(x<0, 0, x)))  # non-negative entries
  P_hat <- t(apply(P2, 1, function(x) x / sum(x)))  # re-scale
  
  # estimation for Theta
  R_hat <- U %*% diag(d) %*% t(V)
  T_hat <- t(solve(t(P_hat) %*% P_hat) %*% t(P_hat) %*% R_hat)
  
  return(list(P_hat=P_hat, T_hat=T_hat, R_hat=R_hat, S_hat=S_hat))
}

