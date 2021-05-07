library(tidyverse)
library(data.table)
library(stringr)
library(future.apply)
library(flock)

#' Master Func function
#'
#' This function Runs the simulation and assigns differet phenotypes (GWAX, LTFH) and generete relevant files used in furhter analysis
#' @param total_indiv The total number of individuals to simulate
#' @param indiv_chunk The number of individuals to simulate in each iteration (should not be above 5000)
#' @param SNP The number of SNPs to simulate
#' @param h the h^2 - meaning the .. ususally 0.5
#' @param c The number of causal SNPs in the simulaiton usually 1/10000*SNP
#' @param k Significance level - usually 0.05
#' @keywords Master
#' @export
#' @examples
#' MasterFunc(10000,1000,10000,0.5,100,0.05)

MasterFunc <- function(total_indiv, indiv_chunk, SNP, h, c, k){
  #Running the simulation
  res <- SimulerData(total_indiv = total_indiv,indiv_chunk = indiv_chunk, SNP = SNP, h = h, c = c, k = k,nr_workers=2)
  #Making a new TRUE file, with correct sibling distribution
  ny_true <- KillSib(fread(paste("./TRUE","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep="")))
  fwrite(as.data.table(ny_true),
         paste("./TRUE","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep=""),
         quote = F,
         sep = " ",
         col.names = T)
  #Making a pheno-file for LTFH with sib
  ltfh_pheno_1 <- Assign_LTFH(Pheno_data = ny_true,
                              valT = k, h2 = h)
  print("LOLL")
  (lfth_pheno_1_selected <- ltfh_pheno_1 %>%
    select(FID,IID,Pheno))

  fwrite(as.data.table(ltfh_pheno_1_selected),
         paste("./Pheno_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1",".txt", sep=""),
         quote = F,
         sep = " ",
         col.names = T)
  #Making a pheno_file for LTFH without sib
  ltfh_pheno_0 <- Assign_LTFH(Pheno_data = ny_true,
                              valT = k, h2 = h, with_sib = 0)
  lfth_pheno_0_selected <- ltfh_pheno_0 %>%
    select(FID,IID,Pheno)

  fwrite(as.data.table(ltfh_pheno_0_selected),
         paste("./Pheno_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0",".txt", sep=""),
         quote = F,
         sep = " ",
         col.names = T)


  #Making histogram-files LTFH with sib
  fwrite(as.data.table(ltfh_pheno_1),
         paste("./Dist_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1",".txt", sep=""),
         quote = F,
         sep = " ",
         col.names = T)

  #Making histogram-files LTFH without sib
  fwrite(as.data.table(ltfh_pheno_0),
         paste("./Dist_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0",".txt", sep=""),
         quote = F,
         sep = " ",
         col.names = T)

  #Making a pheno_file for GWAX with sib
  gwax_pheno_1 <- Assign_GWAX(true = ny_true)
  fwrite(as.data.table(gwax_pheno_1),
         paste("./Pheno_GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1",".txt", sep=""),
         quote = F,
         sep = " ",
         col.names = T)

  #Making a pheno_file for GWAX with out sib
  gwax_pheno_0 <- Assign_GWAX(true = ny_true, with_sib = 0)
  fwrite(as.data.table(gwax_pheno_0),
         paste("./Pheno_GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0",".txt", sep=""),
         quote = F,
         sep = " ",
         col.names = T)

}
