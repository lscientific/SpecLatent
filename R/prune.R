## File Name: prune.R
## File Version: 0.01

#' @title pruning 
#'
#' @description Locate noisy points in the data simplex.
#'
#' @param mat a numeric matrix to be pruned. 
#' @param r the number of neighbors to consider. Default value is 10.
#' @param q the cutoff for the upper quantile of row norms. Higher `q` leads to more points being pruned. Default value is 0.4.
#' @param e the cutoff for the upper quantile of average distance. Higher `e` leads to more points being pruned. Default value is 0.2.
#'
#' @return 
#' \describe{
#'   \item{\code{indices}    the index vector of the rows to be pruned from the left singular matrix.}
#' }
#' 
#' @references Mao, X., Sarkar, P., & Chakrabarti, D. (2021). Estimating mixed memberships with sharp eigenvector deviations. Journal of the American Statistical Association, 116(536), 1928-1940.
#'
#' @export

pruning <- function(mat, r=10, q=0.4, e=0.2) {
  # l2 norms of the rows
  row_norms <- apply(mat, 1, function(x) sqrt(sum(x^2)))
  # indices of points with norm larger than the cutoff
  S_ <- which(row_norms > quantile(row_norms, 1-q))

  x <- c()  # x records the average distance to the r neighbors
  for (s in S_){
    # re-center based on U_{s,:}
    mat_s <- t(apply(mat, 1, function(x) x - mat[s,]))
    # l2 distances to U_{s, :}
    norm_s <- apply(mat_s, 1, function(x) sqrt(sum(x^2)))
    # distance to the r closest neighbors
    d <- norm_s[sort(norm_s, index.return=T)$ix[2:(r+1)]]
    x <- c(x, mean(d))
  }
  indices <- S_[which(x > quantile(x, 1-e))]

  return(indices)
}

