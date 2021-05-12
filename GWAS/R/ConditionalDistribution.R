#' ConditionalDistribution function
#'
#' This function calculates the mean constant and the variance of the normal conditional distributions. Output is a list of the two vector/matrix
#' @param CovMatrix Is the covariance matrix of the multivariate normal distribution
#' @keywords Conditional Distributions
#' @export
#' @examples
#' ConditionalDistribution(matrix(ncol = 4, nrow = 4))

ConditionalDistribution <- function(CovMatrix){
  mu <- matrix(ncol = (nrow(CovMatrix)-1), nrow = nrow(CovMatrix))
  sigma <- numeric(nrow(CovMatrix))
  for (i in 1:nrow(CovMatrix)){
    trans_matrix <- CreateTransMatrix(i,CovMatrix)
    CM <- trans_matrix %*% CovMatrix %*% t(trans_matrix)
    mu[i,] <- 0 + CM[1,2:nrow(CM)] %*% solve(CM[2:nrow(CM),2:nrow(CM)])
    sigma[i] <- CM[1,1] - CM[1,2:nrow(CM)] %*% solve(CM[2:nrow(CM),2:nrow(CM)]) %*% CM[2:nrow(CM),1]
  }
  result <- list("mu" = mu, "sigma" = sigma)
  return(result)
}
