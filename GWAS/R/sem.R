#' Standard Error of the Mean
#'
#' This function calculates the standard error of a vector of means
#' @param x vector of numbers
#' @keywords Standard Error
#' @export
#' @examples
#' sem(rnorm(100))

sem <- function(x){
  sd(x)/sqrt(length(x))
}
