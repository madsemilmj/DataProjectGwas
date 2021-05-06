#' CreateTransMatrix function
#'
#' This function creates a transformation matrix used to calculate conditional distributions
#' @param index Defines the index in the vector we will bring to the top
#' @param covM Defines the CovarianceMatrix for the distribution
#' @keywords transformation
#' @export
#' @examples
#' CreateTransMatrix(1,CovM)

CreateTransMatrix <- function(index, covM){
  trans <- matrix(0, nrow(covM) , nrow(covM))
  diag(trans) <- 1
  res <- trans[-index,]
  res <- rbind(trans[index,],res)
  return(res)
}