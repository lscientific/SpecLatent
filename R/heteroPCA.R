## File Name: heteroPCA.R
## File Version: 0.01


#' @title HeteroPCA 
#' @description This function implements the HeteroPCA algorithm
#'
#' @param R Numeric matrix. The matrix to perform HeteroPCA.
#' @param K Positive integer. The number of top eigenvectors to be extracted.
#' @param T0 Positive integer. The number of iterations.
#' @return Numeric matrix \code{U_hat}
#' @references Zhang, Anru R., T. Tony Cai, and Yihong Wu. "Heteroskedastic PCA: Algorithm, optimality, and applications." The Annals of Statistics 50.1 (2022): 53-80.
#'
#' @export

heteroPCA <- function(R, K, T0) {
  M <- diag_deletion(R %*% t(R))
  M_no_diag <- M
  if (T0 > 0) {
    for (t in 1:T0) {
      svd_res <- svds(M, K)
      if (K == 1) {
        M_bar <- svd_res$u %*% t(svd_res$v) * svd_res$d
      } else {
        M_bar <- svd_res$u %*% diag(svd_res$d[1:K]) %*% t(svd_res$v)
      }
      M <- M_no_diag + diag(diag(M_bar))
    }
  }
  U_hat <- svds(M, K)$u
  return(U_hat)
}
