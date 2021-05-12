#' KillSib function
#'
#' This function assigns the number of siblings to each family-id based on a maksimum likelihood estimate for data on DST.dk
#' @param true Is a TRUE-file from simulation.
#' @keywords Siblings
#' @export
#' @importFrom dplyr %>%
#' @examples
#' KillSib(true = tibble::tibble(ID = 1:100, Child = rbinom(100,1,0.3), Mom = rbinom(100,1,0.3), Dad = rbinom(100,1,0.3), Nr_sib = rep(0,100), Sib_status = rep(0,100), Sib1 = rep(0,100), Sib2 = rep(0,100), Sib3 = rep(0,100), Sib4 = rep(0,100), Sib5 = rep(0,100), Sib6 = rep(0,100), Sib7 = rep(0,100)))


KillSib <- function(true){
  true_new <- true %>%
    tibble::as_tibble() %>%
    dplyr::mutate_all(as.numeric) %>%
    dplyr::mutate(Nr_sib = rbinom(dplyr::n(),7,0.16)) %>%
    tidyr::gather(key,val, Sib1:Sib7) %>%
    dplyr::group_by(ID) %>%
    dplyr::mutate(Sib_status = ifelse(Nr_sib[1]>0,sum(val[c(1:(Nr_sib[1]))]),0)) %>%
    tidyr::spread(key,val)
}
