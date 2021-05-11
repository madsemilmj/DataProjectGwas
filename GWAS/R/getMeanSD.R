library(data.table)
library(tidyverse)
#' Get Mean and SD function
#'
#' This function grabs the mean and sd from the genetic-liability sampled using the Assign_LTFh fuction
#' @param total_indiv Number of indicviduls in the dataset
#' @param SNP Is the number of SNPs in the data
#' @param h Is the heritability usually 0.5
#' @param Child_i Is a binary for indicating case/notcase for child
#' @param Mom_i Is a binary for indicating case/notcase for Mom
#' @param Dad_i Is a binary for indicating case/notcase for Dad
#' @param Nr-Sib_i indiacates the number of siblings for the family
#' @param Sib_status_i indicates the number of sib-cases
#' @keywords Mean/SD
#' @export
#' @examples
#' getMeanSD(total_indiv, SNP, h, Child, Mom, Dad, Nr_sib=0, Sib_status=0)

getMeanSD <- function(total_indiv, SNP, h, Child_i, Mom_i, Dad_i, Nr_sib_i=0, Sib_status_i=0){
  with_sib <- ifelse(Nr_sib_i==0,0,1)
  file_stringer <- paste("./data/Dist_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5_",with_sib,".txt", sep="")
  dist_file <- fread(file_stringer)
  df <- dist_file %>%
    filter(Child == Child_i & Mom == Mom_i & Dad == Dad_i & Nr_sib == Nr_sib_i & Sib_status == Sib_status_i)%>%
    select(Pheno, SDs)
  if (nrow(df)==0){
    #print("FUCK MAN!")
    return(-1)
  }
  else {
    means <- df$Pheno
    SDss <- df$SDs
    mean <- means[1]
    SD <- SDss[1]
    return(list(mean = mean, sd = SD))

  }
}
