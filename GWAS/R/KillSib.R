library(tidyverse)
#' KillSib function
#'
#' This function assigns the number of siblings to each family-id based on a maksimum likelihood estimate for data on DST.dk
#' @param true Is a TRUE-file from simulation.
#' @keywords Siblings
#' @export
#' @examples
#' KillSib(true)


KillSib <- function(true){
  true_new <- true %>%
    as_tibble() %>%
    mutate_all(as.numeric) %>%
    mutate(Nr_sib = rbinom(n(),7,0.16)) %>%
    gather(key,val, Sib1:Sib7) %>%
    group_by(ID) %>%
    mutate(Sib_status = ifelse(Nr_sib[1]>0,sum(val[c(1:(Nr_sib[1]))]),0)) %>%
    spread(key,val)
}
