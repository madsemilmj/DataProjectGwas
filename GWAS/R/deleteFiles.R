#' deleteFiles function
#'
#' This function deletes relevant files
#' @param total_indiv Is the number of individuals for the data-file
#' @param SNP Is the number of SNPs for the data-file
#' @param h Is the heritability usually 0.5
#' @param k siginificance level - set to 0.05
#' @keywords Make Bed
#' @export
#' @examples
#' deleteFiles(total_indiv = 100, SNP = 100, h = 0.5, k = 0.05)


deleteFiles <- function(total_indiv, SNP, h, k){
  files <- character(26)
  files[1] <- paste("./BETA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep="")
  files[2] <- paste("./case_ctrl","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".log", sep="")
  files[3] <- paste("./case_ctrl","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".assoc", sep="")
  files[4] <- paste("./DATA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".bed", sep="")
  files[5] <- paste("./DATA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".bim", sep="")
  files[6] <- paste("./DATA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".fam", sep="")
  files[7] <- paste("./DATA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".log", sep="")
  files[8] <- paste("./DATA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".map", sep="")
  files[9] <- paste("./DATA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".ped", sep="")
  files[10] <- paste("./Dist_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0.txt", sep="")
  files[11] <- paste("./Dist_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1.txt", sep="")
  files[12] <- paste("./GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0.log", sep="")
  files[13] <- paste("./GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0.assoc", sep="")
  files[14] <- paste("./GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1.log", sep="")
  files[15] <- paste("./GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1.assoc", sep="")
  files[16] <- paste("./LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0.log", sep="")
  files[17] <- paste("./LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0.qassoc", sep="")
  files[18] <- paste("./LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1.log", sep="")
  files[19] <- paste("./LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1.qassoc", sep="")
  files[20] <- paste("./MAF","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep="")
  files[21] <- paste("./Pheno","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep="")
  files[22] <- paste("./Pheno_GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0.txt", sep="")
  files[23] <- paste("./Pheno_GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1.txt", sep="")
  files[24] <- paste("./Pheno_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0.txt", sep="")
  files[25] <- paste("./Pheno_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1.txt", sep="")
  files[26] <- paste("./TRUE","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep="")
  files[27] <- paste("./TG","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".log", sep="")
  files[28] <- paste("./TG","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".qassoc", sep="")
  files[29] <- paste("./Pheno_TG","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep="")
  for (file in files){
    if (file.exists(file)){
      file.remove(file)
    }
    else {
      print(paste("Could not find file: ", file, sep=""))
    }
  }
}
