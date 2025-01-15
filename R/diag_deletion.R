## File Name: diag_deletion.R
## File Version: 0.01


#' @title Diagonal deletion 
#' @description This function takes in a matrix, and returns the diagonal-deleted matrix
#'
#' @param X Numeric matrix
#' @return A matrix of with diagonals set to 0
#' @export

diag_deletion <- function(X) {
  if (!("matrix" %in% class(X))) {
    stop("Error: input of diagonal deletion needs to be a matrix!")
  }
  if (nrow(X) != ncol(X)) {
    stop("Error: input of diagonal deletion needs to be a square matrix!")
    
  } 

  diag(X) <- 0
  return(X)
}
