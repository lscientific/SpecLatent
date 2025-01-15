## File Name: DhLCM.R
## File Version: 0.01


#' @title DhLCM clustering and estimation
#' @description This function performs k-means clustering on the top \code{K} eigenvectors/left singular vectors, and estimates the DhLCM model parameters
#'
#' @param R Numeric matrix. Data matrix.
#' @param K Positive integer. The number of top eigenvectors/left singular vectors to be extracted.
#' @param spectral Numeric matrix or character. One of data matrix, \code{"heteroPCA"} and \code{"SVD"}. 
#'        If is a matrix, it is treated as U. Otherwise needs to be a string that 
#'        specifies the method to be used to obtain the top \code{K} 
#'        eigenvectors/left singular vectors. 
#'        \code{"heteroPCA"} implements the heteroPCA method. 
#'        \code{"SVD"} performs ordinary singular vector decomposition.
#' @param norm Character or \code{NULL}. One of \code{"L2"}, \code{"L1"}, \code{"SCORE"}, and \code{NULL}. 
#'        Specifies the method to be used for normalization on the eigenvectors/left singular vectors. 
#'        \code{"L2"} performs L2 normalization. 
#'        \code{"L1"} performs L1 normalization.
#'        \code{"SCORE"} performs SCORE normalization.
#'        \code{"NA"} does not perform normalization.
#' @param dist Character. One of \code{"Bern"}, \code{"Binom"}, and \code{"Pois"}. Specifies the data distribution.
#'        \code{"Bern"} assumes the Bernoulli distribution.
#'        \code{"Binom"} assumes the Binomial distribution.
#'        \code{"Pois"} assumes the Poisson distribution.
#' @param T0 Positive integer. The number of iterations for heteroPCA. Only used when spectral is \code{'heteroPCA'}
#' @param nstart Positive integer. The number of initial starts in the kmeans function.
#' @param S0 Vector or \code{NULL}. If is not \code{NULL}, used to permute the labels.
#' @param clustering_only Boolean. When \code{true}, only clustering is conducted.
#' @return Named list. The list is made of:
#' \itemize{
#' \item \code{U} --- Numeric matrix. Estimation of the left singular matrix.
#' \item \code{T_hat} --- Numeric matrix. Estimation of the \eqn{\Theta}{Theta} matrix.
#' \item \code{sigma2_hat} --- Numeric vector (>=0). Asymptotic variance for each element of \code{T_hat}.
#' \item \code{S_hat} --- Numeric vector. Clustered membership for each subject.
#' \item \code{Z_hat} --- Numeric matrix. Clustered membership for each subject in binary matrix form.
#' }
#' @references Lyu, Zhongyuan, Ling Chen, and Yuqi Gu. "Degree-heterogeneous Latent Class Analysis for High-dimensional Discrete Data." arXiv preprint arXiv:2402.18745 (2024).
#' @export

DhLCM <- function(R, K, spectral='heteroPCA', norm='L2', dist='Bern', 
                  T0=20, nstart=10, S0=NULL, clustering_only=F) { 
  if (!is.matrix(R)) stop("Error: R needs to be an matrix")
  N <- nrow(R)
  J <- ncol(R)
  
  if (is.matrix(spectral)) { # is spectral is an matrix, use it as U
    if ((nrow(spectral) != N) | (ncol(spectral) != K)) {
      stop("Error: spectral matrix dimension incorrect")
    }
    U <- spectral
  } else if (spectral == 'SVD') {
    U <- svds(R, k=K)$u
  } else if (spectral == 'heteroPCA') {
    U <- heteroPCA(R, K, T0)
  } else {
    stop("Error: spectral needs to be an NxK matrix or 'heteroPCA' or 'SVD'")
  }
  
  # k-means clustering
  if (is.null(norm)) {
    kmeans_res <- kmeans(U, centers=K, iter.max=100, nstart=nstart)
  } else if (norm == 'L2') {
    U_bar <- t(apply(U, 1, function(u) u / sqrt(sum(u^2))))
    kmeans_res <- kmeans(U_bar, centers=K, iter.max=100, nstart=nstart)
  } else if (norm == 'L1') {
    U_bar <- t(apply(U, 1, function(u) u / sum(abs(u))))
    kmeans_res <- kmeans(U_bar, centers=K, iter.max=100, nstart=nstart)
  } else if (norm == 'SCORE') {
    U_bar <- t(apply(U, 1, function(u) u / u[1]))
    U_bar <- U_bar[, 2:ncol(U_bar)]
    if (N <= 300) {
      t <- 2 * log(N)
    } else {
      t <- log(N)
    }
    U_bar <- ifelse(abs(U_bar) < t, U_bar, t)
    kmeans_res <- kmeans(U_bar, centers=K, iter.max=100, nstart=nstart)
  } else {
    stop("Error: norm needs to be one of NULL, 'L2', 'L1', 'SCORE'")
  }
  S_hat <- kmeans_res$cluster
  centers <- kmeans_res$centers
  
  # find the optimal permutation if S0 is available
  if (!is.null(S0)) {
    if (length(S0) != N) {
      stop("Error: The length of S0 is not equal to N.")
    }
    
    count1 <- table(S_hat)
    count2 <- table(S0)
    mat <- matrix(0, K, K)
    for (i in 1:K) {
      for (j in 1:K) {
        mat[i, j] <- sum((S_hat == i) & (S0 == j))
      }
    }
    perm_mat <- HungarianSolver(-mat)$pairs[, 2]
    S_hat <- perm(S_hat, perm_mat)
    centers <- centers[sort(perm_mat, index.return=T)$ix, ] 
  }
  
  # estimation
  T_hat = sigma2_hat = Z_hat <- NULL
  if (!clustering_only) {
    Z_hat <- matrix(0, N, K)
    for (i in 1:N) {
      Z_hat[i, S_hat[i]] <- 1
    }
    
    C_hat <- table(S_hat)
    Omega_hat <- apply(U, 1, function(u) sqrt(sum(u^2))) * sapply(S_hat, function(s) sqrt(C_hat[[s]]))
    T_hat <- t(solve(t(Z_hat) %*% Z_hat) %*% t(Z_hat) %*% diag(1 / Omega_hat) %*% R)
    T_hat <- ifelse(T_hat > 1, 1, T_hat)
    T_hat <- ifelse(T_hat < 0, 0, T_hat)
    
    if ((dist == "Bern") | (dist == "Binom")) {
      sigma2_hat <- matrix(NA, J, K)
      for(j in 1:J) {
        for(k in 1:K) {
          w_k <- Omega_hat[S_hat == k]
          sigma2_hat[j, k] <- T_hat[j, k] * sum((1 - w_k * T_hat[j, k]) / w_k) / C_hat[[k]]^2 
        }
      }
    } else if (dist == "Pois") {
      sigma2_hat <- matrix(NA, J, K)
      for(j in 1:J) {
        for(k in 1:K) {
          w_k <- Omega_hat[S_hat == k]
          sigma2_hat[j, k] <- T_hat[j, k] * sum(1 / w_k)/ C_hat[[k]]^2
        }
      }
    } else {
      stop("Error: dist needs to be 'Bern' or 'Binom' or 'Pois'")
    }
  }
  
  return(list(U=U, T_hat=T_hat, sigma2_hat=sigma2_hat, S_hat=S_hat, Z_hat=Z_hat))
}


