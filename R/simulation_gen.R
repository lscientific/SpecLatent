## File Name: simulation_gen.R
## File Version: 0.01
## generates simulated data

library(gtools)
library(MASS)
library(stats)


#' @title Generates polytomous data 
#' @description This function generates polytomous GoM data
#' @param N integer. The number of subjects
#' @param J integer. The number of items
#' @param K integer. The number of extreme latent profiles
#' @param C integer. The number of possible outcomes for each item
#' @param alphas_pi vector of numeric values. Specifies the Dirichlet parameters for membership score generation. Need to be length of \code{K}
#' #' @param alphas_theta vector of numeric values. Specifies the Dirichlet parameters for item parameter generation. Need to be length of \code{C}
#' @return Named list. The list is made of:
#' \itemize{
#' \item \code{Pi} --- Numeric matrix. Membership score matrix
#' \item \code{Theta} --- Numeric matrix. Flattened item parameter matrix
#' \item \code{R_flattened} --- Binary matrix of size \code{N}x(\code{CJ}). Generated flattened data matrix 
#' \item \code{R_pol} --- Integer matrix of size \code{N}x\code{J}. Generated polytomous data matrix 
#' }
#' @export
#' 

data_gen_pol <- function(N, J, K, C, alphas_pi, alphas_theta) {
  if (length(alphas_pi) != K) stop("Error: the length of alphas_pi needs to K!")
  if (length(alphas_theta) != C) stop("Error: the length of alpha_theta needs to C!")
  
  # membership score matrix
  Pi <- rdirichlet(N, alphas_pi)
  Pi[1:K, ] <- diag(1, K)
  
  # flattened item parameter matrix
  Theta <- matrix(NA, nc=K, nr=J*C)
  for(j in 1:J) {
    for (k in 1:K) {
      Theta[((j-1)*C+1):(j*C), k] <- rdirichlet(1, alphas_theta)
    }
  }
  
  # flattened data matrix
  R_flattened <- matrix(NA, N, C*J)
  R_pol <- matrix(NA, N, J)
  for(i in 1:N) {
    for(j in 1:J) {
      prob <- matrix(Pi[i, ], nr=1) %*% t(Theta[((j-1)*C+1):(j*C), ])
      response <- rmultinom(1, 1, prob)
      R_flattened[i, ((j-1)*C+1):(j*C)] <- response
      R_pol[i, j] <- which(response == 1) - 1
    }
  }
  
  return(list(Pi=Pi, Theta=Theta, R_flattened=R_flattened, R_pol=R_pol))
}



#' @title Generates locally dependent Bernoulli data 
#' @description This function generates locally dependent Bernoulli distributed GoM data.
#' @param N integer. The number of subjects
#' @param J integer. The number of items
#' @param K integer. The number of extreme latent profiles
#' @param rho numeric. The auto-regressive parameter for local dependence. Needs to be in \code{[0,1]}
#' @param block integer. The local dependence block size
#' @param alphas_pi vector of numeric values. Specifies the Dirichlet parameters for membership score generation. Need to be length of \code{K}
#' @param pars vector of non-negative numeric values. Specifies the Beta distribution shape parameters for item parameter generation. Need to be length of \code{2}
#' @return Named list. The list is made of:
#' \itemize{
#' \item \code{Pi} --- Numeric matrix. Membership score matrix
#' \item \code{Theta} --- Numeric matrix. Flattened item parameter matrix
#' \item \code{R} --- Binary matrix of size \code{N}x\code{J}. Generated locally dependent data matrix 
#' }
#' @export

data_gen_loc <- function(N, J, K, rho, block, alphas_pi, pars=c(0.2, 0.2)) {
  if (length(alphas_pi) != K) stop("Error: the length of alphas_pi needs to K!")
  if (length(pars) != 2) stop("Error: the length of pars needs to 2!")
  if (any(pars < 0)) stop("Error: pars needs to be non-negative!")
  
  # membership score matrix
  Pi <- rdirichlet(N, alphas_pi)
  Pi[1:K, ] <- diag(1, K)
  
  # item parameter matrix
  Theta <- matrix(rbeta(J*K, pars[1], pars[2]), nr=J, nc=K)
  
  # Normal quantile
  R0 <- Pi %*% t(Theta)
  D <- matrix(NA, nrow=N, ncol=J)
  for(j in 1:J) {
    D[, j] <- sapply(R0[, j], qnorm)
  }
  
  # correlation matrix
  exponent <- abs(matrix(1:block-1, nrow=block, ncol=block, byrow=TRUE) - (1:block-1))
  Sigma <- rho^exponent
  
  # data matrix
  R <- matrix(NA, nrow=N, ncol=J)
  for(i in 1:N) {
    Y <- mvrnorm(J/block, mu=rep(0, block), Sigma=Sigma)
    Y <- as.vector(t(Y))
    R[i, ] <- as.numeric(Y < D[i, ])
  }
  
  return(list(Pi=Pi, Theta=Theta, R=R))
}




#' @title Generates Poisson data 
#' @description This function generates Poisson distributed GoM data.
#' @param N integer. The number of subjects
#' @param J integer. The number of items
#' @param K integer. The number of extreme latent profiles
#' @param alphas_pi vector of numeric values. Specifies the Dirichlet parameters for membership score generation. Need to be length of \code{K}
#' @param pars vector of positive numeric values. Specifies the Gamma distribution shape and rate parameters for item parameter generation. Need to be length of \code{2}
#' @return Named list. The list is made of:
#' \itemize{
#' \item \code{Pi} --- Numeric matrix. Membership score matrix
#' \item \code{Theta} --- Numeric matrix. Flattened item parameter matrix
#' \item \code{R} --- Binary matrix of size \code{N}x\code{J}. Generated Poisson data matrix 
#' }
#' @export

data_gen_pois <- function(N, J, K, alphas_pi, pars=c(1, 2)) {
  if (length(alphas_pi) != K) stop("Error: the length of alphas_pi needs to K!")
  if (length(pars) != 2) stop("Error: the length of pars needs to 2!")
  if (any(pars <= 0)) stop("Error: pars needs to be positive!")
  
  # membership score matrix
  Pi <- rdirichlet(N, alphas_pi)
  Pi[1:K, ] <- diag(1, K)
  
  # item parameter matrix
  Theta <- matrix(rgamma(J*K, pars[1], pars[2]), nr=J, nc=K)
  
  # data matrix
  R0 <- Pi %*% t(Theta)
  R <- matrix(NA, nrow=N, ncol=J)
  for(i in 1:N) {
    for(j in 1:J) {
      R[i, j] <- rpois(1, R0[i,j])
    }
  }
  
  return(list(Pi=Pi, Theta=Theta, R=R))
}


