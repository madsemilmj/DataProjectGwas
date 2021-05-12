#' covariance matrix function
#'
#' This function calculates the covariance matrix based on h^2 and the number of siblings
#' @param h2 is the .. usually 0.5
#' @param nr_sib Defines the number of sibling, defaults to 0
#' @keywords covariance matrix
#' @export
#' @examples
#' CovarianceMatrix(0.5)

CovarianceMatrix <- function(h2, nr_sib = 0) {
  cov <- matrix(h2/2, 4+nr_sib , 4+nr_sib)
  diag(cov) <- 1
  cov[2:(4+nr_sib),1] <- cov[1,2:(4+nr_sib)] <- cov[4, 3] <- cov[3,4] <-0
  cov[2, 2] <- h2
  cov[1, 1] <- 1-h2
  cov
}
