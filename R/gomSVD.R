## File Name: gomSVD.R
## File Version: 0.01


#' @title gomSVD
#'
#' @description Estimation algorithm for generalized-GoM model with locally dependent data
#'
#' @param U the pruned left singular matrix from data SVD. 
#' @param V the right singular matrix from data SVD. 
#' @param d the vector containing the singular values. 
#' @param prune logical; if true, the pruning step is performed.
#' @param r the number of neighbors to consider in pruning. Used only when `prune=T`. Default value is 10.
#' @param q the cutoff for the upper quantile of row norms. Used only when `prune=T`. Higher `q` leads to more points being pruned. Default value is 0.4.
#' @param e the cutoff for the upper quantile of average distance. Used only when `prune=T`. Higher `e` leads to more points being pruned. Default value is 0.2.
#' @param lower the minimum value for item parameters. Default value is 0.
#' @param upper the minimum value for item parameters. Default value is 1.
#' @return The function returns a list with the following components:
#' \describe{
#'   \item{`P_hat`}{The estimated membership scores.}
#'   \item{`T_hat`}{The estimated item response parameters.}
#'   \item{`R_hat`}{The estimated response expectation.}
#'   \item{`S_hat`}{The estimated indices of pure subjects.}
#'   \item{`t`}{Computation time.}
#' }

gomSVD <- function(U, V, d, lower=0, upper=1,
                    prune=T, r=10, q=0.4, e=0.2) {
  t1 <- Sys.time()
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
  
  # threshold
  R_hat[R_hat < lower] <- lower
  R_hat[R_hat > upper] <- upper
  T_hat[T_hat < lower] <- lower
  T_hat[T_hat > upper] <- upper
  
  t2 <- Sys.time()
  
  return(list(P_hat=P_hat, T_hat=T_hat, R_hat=R_hat, S_hat=S_hat, t=t2-t1))
}

