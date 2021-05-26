#' GwasPlink Function for mac
#'
#' This function converts the ped-file to the less of size bed file (used by PLINK)
#' @param total_indiv Is the number of individuals for the data-file
#' @param SNP Is the number of SNPs for the data-file
#' @param h Is the heritability usually 0.5
#' @param k siginificance level - set to 0.05
#' @keywords Make Bed
#' @export
#' @examples
#' gwasPlinkmac(total_indiv = 1000, SNP = 1000, h = 0.5, k = 0.05)


gwasPlinkmac <- function(total_indiv, SNP, h, k){
  bfile <-paste("DATA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100, sep="")
  phenotype <- paste("Pheno","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep="")
  phenotype_GWAX_0 <- paste("Pheno_GWAX_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0",".txt", sep="")
  phenotype_GWAX_1 <- paste("Pheno_GWAX_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1",".txt", sep="")
  phenotype_LTFH_0 <- paste("Pheno_LTFH_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0",".txt", sep="")
  phenotype_LTFH_1 <- paste("Pheno_LTFH_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1",".txt", sep="")
  phenotype_TG <- paste("Pheno_TG_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep="")
  case_ctrl <- paste("case_ctrl","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100, sep="")
  GWAX_0 <- paste("GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0", sep="")
  GWAX_1 <- paste("GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1", sep="")
  LTFH_0 <- paste("LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0", sep="")
  LTFH_1 <- paste("LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1", sep="")
  TRUE_GEN <- paste("TG","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100, sep="")
  makegwas1<- paste("plink --bfile", bfile, "--pheno", phenotype_GWAX_0, "--pheno-name Pheno --out", GWAX_0, "--assoc")
  makegwas2<- paste("plink --bfile", bfile, "--pheno", phenotype_GWAX_1, "--pheno-name Pheno --out", GWAX_1, "--assoc")
  makegwas3<- paste("plink --bfile", bfile, "--pheno", phenotype_LTFH_0, "--pheno-name Pheno --out", LTFH_0, "--assoc")
  makegwas4<- paste("plink --bfile", bfile, "--pheno", phenotype_LTFH_1, "--pheno-name Pheno --out", LTFH_1, "--assoc")
  makegwas5<- paste("plink --bfile", bfile, "--pheno", phenotype, "--pheno-name Pheno --out", case_ctrl, "--assoc")
  makegwas6<- paste("plink --bfile", bfile, "--pheno", phenotype_TG, "--pheno-name Pheno --out", TRUE_GEN, "--assoc")
  system(command = makegwas1)
  system(command = makegwas2)
  system(command = makegwas3)
  system(command = makegwas4)
  system(command = makegwas5)
  system(command = makegwas6)
}
