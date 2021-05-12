#' Assing GWAX function
#'
#' This function assigns phenotypes based on GWAX.
#' @param true Is a TRUE-file from simulation.
#' @param with_sib Defines if siblings Default is 1=true, else 0
#' @keywords gwax
#' @export
#' @importFrom dplyr %>%
#' @examples
#' Assign_GWAX(tibble::tibble(ID = 1:100, Child = rbinom(100,1,0.3), Mom = rbinom(100,1,0.3), Dad = rbinom(100,1,0.3), Nr_sib = rep(0,100), Sib_status = rep(0,100)), with_sib = 1)
Assign_GWAX <- function(true, with_sib = 1){
  res <- true %>%
    tibble::as_tibble()%>%
    dplyr::select(ID, Child, Mom, Dad, Nr_sib, Sib_status)
  if (with_sib == 0){
    res$Nr_sib <- 0
    res$Sib_status <- 0
  }
  res$Pheno <- ifelse(res$Child > 0 | res$Mom > 0 | res$Dad > 0 | res$Sib_status > 0 ,1,0)
  res$FID <- res$ID
  result <- res %>%
    dplyr::rename(IID = ID) %>%
    dplyr::select(FID,IID,Pheno)
  return(result)

}
