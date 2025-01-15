## File Name: gGoM.R
## File Version: 0.01


#' @title Generalized-GoM estimation
#'
#' @description Estimation algorithm for generalized-GoM models
#'
#' @param R data matrix.
#' @param K integer. The number of extreme latent profiles. \code{K} should be at least 2.
#' @param pol logical; if true, assume GoM model with polytomous response, and flattening is applied. Item parameter estimation \code{T_hat} is also flattened.
#' @param large logical; if true, \code{K} needs to be at least 3 and use the large-scale SVD function \code{RSpectra::svds}.
#' @param prune logical; if true, the pruning step is performed.
#' @param dist character; One of \code{"Bern"}, \code{"Binom"}, and \code{"Pois"}. Specifies the data distribution.
#'        \code{"Bern"} assumes the Bernoulli distribution.
#'        \code{"Binom"} assumes the Binomial distribution.
#'        \code{"Pois"} assumes the Poisson distribution.
#' @param r the number of neighbors to consider in pruning. Used only when \code{prune} is \code{TRUE}. Default value is 10.
#' @param q the cutoff for the upper quantile of row norms. Used only when \code{prune} is \code{TRUE}. Higher \code{q} leads to more points being pruned. Default value is 0.4.
#' @param e the cutoff for the upper quantile of average distance. Used only when \code{prune} is \code{TRUE}. Higher \code{e} leads to more points being pruned. Default value is 0.2.
#' @param lower the minimum value for item parameters. Default value is 0.
#' @param upper the minimum value for item parameters. Default value is 1.
#' @return The function returns a list with the following components:
#' \itemize{
#' \item{\code{P_hat} the estimated membership scores.}
#' \item{\code{T_hat} the estimated item response parameters.}
#' \item{\code{R_hat} the estimated response expectation.}
#' \item{\code{S_hat} the estimated indices of pure subjects.}
#' \item{\code{t} computation time.}
#' }
#' @references Chen, Ling, and Yuqi Gu. "A spectral method for identifiable grade of membership analysis with binary responses." Psychometrika (2024): 1-32.
#' @references Chen, Ling, Chengzhu Huang, and Yuqi Gu. "Generalized Grade-of-Membership Estimation for High-dimensional Locally Dependent Data." arXiv preprint arXiv:2412.19796 (2024).
#' @export

gGoM <- function(R, K, pol=T, dist=NULL, large=T, prune=T, r=10, q=0.4, e=0.2) {
  t1 <- Sys.time()
  
  if (pol) {
    # the number of categories for each item
    Cs <- sapply(apply(R, 2, function(x) sort(unique(x))), length)
    # flatten the polytomous matrix into a `fat` binary matrix
    R <- flatten(R)
  } else {
    if (is.null(dist)) {
      stop("Error: `dist` input required when `pol` is FALSE.")
    } 
    if ((dist == "Bern") | (dist == "Binom") ) {
      lower <- 0
      upper <- 1
    } else if (dist == "Pois") {
      lower <- 0
      upper <- Inf
    }
  }
  
  # SVD
  if (large) {
    if (K < 3) stop("Error: K should be at least 3 when `large` is TRUE")
    svd_res <- RSpectra::svds(R, k=K)
  } else {
    svd_res <- svd(R, nu=K, nv=K)
  }
  
  # parameter estimation
  res <- gomSVD(svd_res$u, svd_res$v, svd_res$d, prune, r, q, e)
  
  # post-processing
  if (pol) {
    # re-scale the item parameter matrix estimation
    res$T_hat <- rescale_T(res$T_hat, Cs)
  } else {
    if (dist == "Bern") res$T_hat <- res$T_hat / 2
    # threshold
    res$R_hat[res$R_hat < lower] <- lower
    res$R_hat[res$R_hat > upper] <- upper
    T_hat[T_hat < lower] <- lower
    T_hat[T_hat > upper] <- upper
  }
  
  t2 <- Sys.time()
  res$t <- t2-t1
  
  return(res)
}


