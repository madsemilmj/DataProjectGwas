#' MakeBed function for mac
#'
#' This function converts the ped-file to the less of size bed file (used by PLINK)
#' @param total_indiv Is the number of individuals for the data-file
#' @param SNP Is the number of SNPs for the data-file
#' @param h Is the heritability usually 0.5
#' @param k siginificance level - set to 0.05
#' @keywords Make Bed
#' @export
#' @examples
#' makeBedmac(total_indiv = 1000, SNP = 1000, h = 0.5, k = 0.05)

makeBedmac <- function(total_indiv, SNP, h, k){
  pfile <-paste("DATA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100, sep="")
  bfile <-paste("DATA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100, sep="")
  bed <- paste("plink --file", pfile, "--make-bed --out", bfile)
  system(command=bed)
}
