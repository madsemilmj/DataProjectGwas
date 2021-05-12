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
#' @importFrom dplyr %>%
#' @examples
#' MasterFunc(1000,100,1000,0.5,100,0.05)

MasterFunc <- function(total_indiv, indiv_chunk, SNP, h, c, k){
  # Check if data folder is created, if not break simulation
  mainDir <- getwd()
  subDir <- "data"
  does_folder_exists <- ifelse(dir.exists(file.path(mainDir, subDir)),1,0)
  if (does_folder_exists == 1){
    #Running the simulation
    res <- SimulerData(total_indiv = total_indiv,indiv_chunk = indiv_chunk, SNP = SNP, h = h, c = c, k = k,nr_workers=2)
    #Making a new TRUE file, with correct sibling distribution
    ny_true <- KillSib(data.table::fread(paste("./TRUE","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep="")))
    data.table::fwrite(data.table::as.data.table(ny_true),
                       paste("./TRUE","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep=""),
                       quote = F,
                       sep = " ",
                       col.names = T)
    #Making a pheno-file for LTFH with sib
    ltfh_pheno_1 <- Assign_LTFH(Pheno_data = ny_true,
                                valT = k, h2 = h)
    lfth_pheno_1_selected <- ltfh_pheno_1 %>%
      dplyr::select(FID,IID,Pheno)


    data.table::fwrite(data.table::as.data.table(lfth_pheno_1_selected),
                       paste("./Pheno_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1",".txt", sep=""),
                       quote = F,
                       sep = " ",
                       col.names = T)
    #Making a pheno_file for LTFH without sib
    ltfh_pheno_0 <- Assign_LTFH(Pheno_data = ny_true,
                                valT = k, h2 = h, with_sib = 0)
    lfth_pheno_0_selected <- ltfh_pheno_0 %>%
      dplyr::select(FID,IID,Pheno)

    data.table::fwrite(data.table::as.data.table(lfth_pheno_0_selected),
                       paste("./Pheno_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0",".txt", sep=""),
                       quote = F,
                       sep = " ",
                       col.names = T)


    #Making histogram-files LTFH with sib
    data.table::fwrite(data.table::as.data.table(ltfh_pheno_1),
                       paste("./Dist_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1",".txt", sep=""),
                       quote = F,
                       sep = " ",
                       col.names = T)

    #Making histogram-files LTFH without sib
    data.table::fwrite(data.table::as.data.table(ltfh_pheno_0),
                       paste("./Dist_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0",".txt", sep=""),
                       quote = F,
                       sep = " ",
                       col.names = T)

    #Making a pheno_file for GWAX with sib
    gwax_pheno_1 <- Assign_GWAX(true = ny_true)
    data.table::fwrite(data.table::as.data.table(gwax_pheno_1),
                       paste("./Pheno_GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1",".txt", sep=""),
                       quote = F,
                       sep = " ",
                       col.names = T)

    #Making a pheno_file for GWAX with out sib
    gwax_pheno_0 <- Assign_GWAX(true = ny_true, with_sib = 0)
    data.table::fwrite(data.table::as.data.table(gwax_pheno_0),
                       paste("./Pheno_GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0",".txt", sep=""),
                       quote = F,
                       sep = " ",
                       col.names = T)
    #Makebed function
    makeBed(total_indiv, SNP, h, k)

    #Rungwas
    gwasPlink(total_indiv, SNP, h, k)

    #Moverelevantfiles
    moveRelevantFiles(total_indiv, SNP, h, k)

    #Delete files
    deleteFiles(total_indiv, SNP, h, k)

  } else {
    print("Simulation failed - please create a folder named 'data' in your working directory")
  }

}
