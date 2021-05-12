#' moveRelevantFiles function
#'
#' This function moves relevant files
#' @param total_indiv Is the number of individuals for the data-file
#' @param SNP Is the number of SNPs for the data-file
#' @param h Is the heritability usually 0.5
#' @param k siginificance level - set to 0.05
#' @keywords Make Bed
#' @export
#' @examples
#' moveRelevantFiles(total_indiv = 1000, SNP = 1000, h = 0.5, k = 0.05)


moveRelevantFiles <- function(total_indiv, SNP, h, k){
  #Checking if folder exists otherwise abort
  mainDir <- getwd()
  subDir <- "data"
  does_folder_exists <- ifelse(dir.exists(file.path(mainDir, subDir)),1,0)
  if (does_folder_exists == 1){
    #Checking if LTFH qassoc exists since we then assume everything is ok
    if (file.exists(paste("./LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1.qassoc", sep=""))){
      #Moving BETA
      beta <- data.table::fread(paste("./BETA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep=""))
      data.table::fwrite(data.table::as.data.table(beta),
                         paste("./data/BETA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep=""),
                         quote = F,
                         sep = " ",
                         col.names = T)
      #Moving case_control.assoc
      case_control <- data.table::fread(paste("./case_ctrl","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".assoc", sep=""))
      data.table::fwrite(data.table::as.data.table(case_control),
                         paste("./data/case_ctrl","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".assoc", sep=""),
                         quote = F,
                         sep = " ",
                         col.names = T)
      #Moving DIST_LTFH_1
      DIST_LTFH_1 <- data.table::fread(paste("./Dist_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1.txt", sep=""))
      data.table::fwrite(data.table::as.data.table(DIST_LTFH_1),
                         paste("./data/Dist_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1.txt", sep=""),
                         quote = F,
                         sep = " ",
                         col.names = T)

      #Moving DIST_LTFH_0
      DIST_LTFH_0 <- data.table::fread(paste("./Dist_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0.txt", sep=""))
      data.table::fwrite(data.table::as.data.table(DIST_LTFH_1),
                         paste("./data/Dist_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0.txt", sep=""),
                         quote = F,
                         sep = " ",
                         col.names = T)

      #Moving GWAX_0.assoc
      GWAX_0_assoc <- data.table::fread(paste("./GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0.assoc", sep=""))
      data.table::fwrite(data.table::as.data.table(GWAX_0_assoc),
                         paste("./data/GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0.assoc", sep=""),
                         quote = F,
                         sep = " ",
                         col.names = T)

      #Moving GWAX_1.assoc
      GWAX_1_assoc <- data.table::fread(paste("./GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1.assoc", sep=""))
      data.table::fwrite(data.table::as.data.table(GWAX_1_assoc),
                         paste("./data/GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1.assoc", sep=""),
                         quote = F,
                         sep = " ",
                         col.names = T)

      #Moving LTFH_0.qassoc
      LTFH_0_qassoc <- data.table::fread(paste("./LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0.qassoc", sep=""))
      data.table::fwrite(data.table::as.data.table(LTFH_0_qassoc),
                         paste("./data/LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0.qassoc", sep=""),
                         quote = F,
                         sep = " ",
                         col.names = T)

      #Moving LTFH_1.qassoc
      LTFH_1_qassoc <- data.table::fread(paste("./LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1.qassoc", sep=""))
      data.table::fwrite(data.table::as.data.table(LTFH_1_qassoc),
                         paste("./data/LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1.qassoc", sep=""),
                         quote = F,
                         sep = " ",
                         col.names = T)

      #Moving MAF
      MAF <- data.table::fread(paste("./MAF","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep=""))
      data.table::fwrite(data.table::as.data.table(MAF),
                         paste("./data/MAF","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep=""),
                         quote = F,
                         sep = " ",
                         col.names = T)
      #Moving TRUE-file
      true <- data.table::fread(paste("./TRUE","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep=""))
      data.table::fwrite(data.table::as.data.table(true),
                         paste("./data/TRUE","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep=""),
                         quote = F,
                         sep = " ",
                         col.names = T)
    } else {print("file not exists!")}
  }else {
    print("Simulation failed - please create a folder named 'data' in your working directory")
  }

  }
