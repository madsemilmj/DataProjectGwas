library(tidyverse)
#' Assing GWAX function
#'
#' This function assigns phenotypes based on GWAX.
#' @param true Is a TRUE-file from simulation.
#' @param with_sib Defines if siblings Default is 1=true, else 0
#' @keywords gwax
#' @export
#' @examples
#' Assign_GWAX()
Assign_GWAX <- function(true, with_sib = 1){
  res <- true %>%
    as_tibble()%>%
    select(ID, Child, Mom, Dad, Nr_sib, Sib_status)
  if (with_sib == 0){
    res$Nr_sib <- 0
    res$Sib_status <- 0
  }
  res$Pheno <- ifelse(res$Child > 0 | res$Mom > 0 | res$Dad > 0 | res$Sib_status > 0 ,1,0)
  res$FID <- res$ID
  result <- res %>%
    rename(IID = ID) %>%
    select(FID,IID,Pheno)
  return(result)

}

