#' LinearTransformation function
#'
#' This function makes a linear transformation of a CovarianceMatrix.
#' @param covM Is a covariance matrix
#' @keywords Linaer Transformation
#' @export
#' @examples
#' LinearTransformation(covM)
LinearTransformation = function(covM){ 
  trans <- matrix(0, nrow(covM) , nrow(covM) )
  diag(trans) <- 1
  trans[2,1]<- trans[1,2] <-1
  trans[1,1]<- 0
  trans %*% covM %*% t(trans) 
}