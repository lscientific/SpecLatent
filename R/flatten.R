## File Name: flatten.R
## File Version: 0.01


#' @title flatten 
#' @description Flatten the polytomous matrix to a fat binary matrix
#'
#' @param R integer matrix. The polytomous response data matrix. 
#'
#' @return 
#' \describe{
#'   \item{\code{R_flattened}    flattened binary matrix.}
#' }
#' 


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
