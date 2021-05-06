#' Transform to PLINK function
#'
#' This function takes a matrix of 0,1,2 (genotypes) and transforms it to data readable to PLINK-software.
#' @param df Is a matrix/data-frame that we want to Normalize (number_of_indivs X SNPs)
#' @param indiv_chunk Defines the number of individuals in the dataframe
#' @keywords PLINK
#' @export
#' @examples
#' TransformToPlink(df, 1000)

TransformToPlink <- function(df,indiv_chunk){
  a <- matrix(ncol=3,nrow=2) #C er den riskogivne dvs CC = 2 og AA = 0
  a[1,]<-c("A", "C", "C")
  a[2,]<- c("A", "A", "C")
  
  for_plink <- lapply(1:indiv_chunk, function(x){
    tmp <- a[,df[x,]+1]
    c(rep(1,4),rbind(tmp[1,],tmp[2,]))
  })%>%
    do.call(rbind, .)
  return(for_plink)
}