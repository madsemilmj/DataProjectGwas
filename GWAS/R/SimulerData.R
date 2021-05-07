library(tidyverse)
library(data.table)
library(stringr)
library(future.apply)
library(flock)

#' SimulerData function
#'
#' This function simulates genetic-data and saves relevant files doing so.
#' @param total_indiv The total number of individuals to simulate
#' @param indiv_chunk The number of individuals to simulate in each iteration (should not be above 5000)
#' @param SNP The number of SNPs to simulate
#' @param h the h^2 - meaning the .. ususally 0.5
#' @param c The number of causal SNPs in the simulaiton usually 1/10000*SNP
#' @param k Significance level - usually 0.05
#' @param nr_workers Number of cores to run the code in parallel
#' @keywords Genetic simulation
#' @export
#' @examples
#' SimulerData(10000,1000,10000,0.5,100,0.05,2)


SimulerData <- function(total_indiv, indiv_chunk, SNP, h, c, k, nr_workers) {
  #loader <- image_read("./giphy.gif")
  print("Simulation started, DONE! will be printed when the simulaiton is done")
  #print(loader)
  t <- qnorm(1-k, 0, 1)
  MAF <- runif(n = SNP, min = 0.01, max = 0.49) # Simulating #SNP probabilities.
  #Generating MAP file for plink
  MAP <- cbind(1, 1:SNP, 0, 1:SNP)
  fwrite(as.data.table(MAP),
         paste("./DATA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".map", sep=""),
         quote = F,
         sep = " ",
         col.names = F)

  #Writing the probabilities to a txt-file.
  fwrite(as.data.table(matrix(MAF,ncol = length(MAF))),
         paste("./MAF","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep=""),
         quote = F,
         sep = " ",
         col.names = F)

  # Drawing C indices from 1:SNP without replacement (represents the C causal SNPS)
  disease <- sample(SNP, c, replace = FALSE, prob = NULL)

  # Using the above indices to find the effect sizes
  beta <- numeric(SNP)
  beta[disease] <- rnorm(c, 0, sqrt(h/c))


  # Defining genotypes for mom, dad and child
  df_mom <- t(replicate(indiv_chunk, rbinom(SNP, 2, MAF)))
  df_dad <- t(replicate(indiv_chunk, rbinom(SNP, 2, MAF)))
  df_child <- lapply(1:indiv_chunk, function(x){round(((df_mom[x,] + df_dad[x,])/2)+rnorm(SNP,mean = 0, sd = 0.005 ))}) %>%
    do.call(rbind, .)
  sib_matrix <- matrix(nrow= indiv_chunk, ncol = 7)
  for (i in 1:7){
    df_sib <- lapply(1:indiv_chunk, function(x){round(((df_mom[x,] + df_dad[x,])/2)+rnorm(SNP,mean = 0, sd = 0.005 ))}) %>%
      do.call(rbind, .)
    df_sib <- Normalize_data(df_sib, MAF)
    l_g_sib <- df_sib %*% beta
    l_e_sib <- rnorm(indiv_chunk, 0, sqrt(1-h))
    l_sib <- l_g_sib + l_e_sib
    sib_matrix[,i] <- ifelse(l_sib >= t, 1, 0)
  }


  #NORMALIZING MATRICES
  df_mom <- Normalize_data(df_mom, MAF) #We can overwrite because df_mom no longer needed
  df_dad <- Normalize_data(df_dad, MAF) #We can overwrite because df_dad no longer needed
  df_child_normal <- Normalize_data(df_child, MAF)


  ### DEFINING LIABILITIES
  l_g_child <- df_child_normal %*% beta

  #Sampling the environmental liability (child)
  l_e_child <- rnorm(indiv_chunk, 0, sqrt(1-h))

  # Finding the full-liability (child)
  l_child <- l_g_child + l_e_child
  # Using the full-liability to determine binary phenotype (child)
  y_child <- ifelse(l_child >= t, 1, 0)

  ######DAD#########
  # GENETISK LIABILITY (dad) (df_dad is already normalized (overwritten))
  l_g_dad <- df_dad %*% beta
  l_e_dad <- rnorm(indiv_chunk, 0, sqrt(1-h))
  l_dad <- l_g_dad + l_e_dad
  y_dad <- ifelse(l_dad >= t, 1, 0)

  ###### MOM ################

  # GENETISK LIABILITY (mom) (df_mom is already normalized (overwritten))
  l_g_mom <- df_mom %*% beta
  l_e_mom <- rnorm(indiv_chunk, 0, sqrt(1-h))
  l_mom <- l_g_mom + l_e_mom
  y_mom <- ifelse(l_mom >= t, 1, 0)



  true_vals <- cbind(1:indiv_chunk,l_g_child,l_e_child, l_child, y_child, l_g_dad, l_e_dad, l_dad, y_dad, l_g_mom, l_e_mom, l_mom, y_mom, sib_matrix)
  colnames(true_vals) <- c("ID", "LG_child", "LE_child", "L_child","Child","LG_dad", "LE_dad", "L_dad", "Dad", "LG_mom", "LE_mom", "L_mom", "Mom", "Sib1", "Sib2", "Sib3", "Sib4", "Sib5", "Sib6", "Sib7")

  # Writing true vals to file
  fwrite(as.data.table(true_vals),
         paste("./TRUE","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep=""),
         quote = F,
         sep = " ",
         col.names = T)

  # Writing effectsizes to file
  fwrite(as.data.table(matrix(beta,ncol = length(beta))),
         paste("./BETA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep=""),
         quote = F,
         sep = " ",
         col.names = F)

  #Making space in memory
  rm(df_dad,df_mom,df_child_normal,MAP, sib_matrix, df_sib)
  #Filling data in the for_plink matrix (for the first #indiv_chunk persons)
  for_plink <- TransformToPlink(df_child,indiv_chunk)
  for_plink <- cbind(1:indiv_chunk,1:indiv_chunk,for_plink) #id

  #Writing PED file for plink
  fwrite(as.data.table(for_plink),
         paste("./DATA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".ped", sep=""),
         quote = F,
         sep = " ",
         col.names = F)

  # Creating and writing PHENO-types for child
  PHENO <- cbind(1:indiv_chunk,1:indiv_chunk,y_child)
  colnames(PHENO) <- c("FID","IID", "Pheno" )
  fwrite(as.data.table(PHENO),
         paste("./Pheno","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep=""),
         quote = F,
         sep = " ",
         col.names = T)

  rm(for_plink,PHENO)
  #Locking a tmp-file to make sure we are able to write to a file simultaneous when using parallel computing
  lock = tempfile()

  #### you generally do not want to put plan inside a function like this. future_lapply collapses to a normal lapply function if no plan is given.
  #### this allows the end user to use it however they want.
  # Defining the number of cores in the processor to run simultaneous
  x <<- nr_workers #GLOBAL VARIABLE SO IT CAN BE USED BY FUTURE
  plan(future::multisession(workers = x))
  nr_it <- total_indiv %/% indiv_chunk #Integer division
  #Performing the steps in parallel
  tmp<- future.apply::future_lapply(2:nr_it, function(n) {

    id_fra <- (indiv_chunk*(n-1))+1
    id_til <- ((n-1)+1)*indiv_chunk

    # Defining genotypes for mom, dad and child
    df_mom <- t(replicate(indiv_chunk, rbinom(SNP, 2, MAF)))
    df_dad <- t(replicate(indiv_chunk, rbinom(SNP, 2, MAF)))
    df_child <- lapply(1:indiv_chunk, function(x){round(((df_mom[x,] + df_dad[x,])/2)+rnorm(SNP,mean = 0, sd = 0.005 ))}) %>%
      do.call(rbind, .)
    sib_matrix <- matrix(nrow= indiv_chunk, ncol = 7)
    for (i in 1:7){
      df_sib <- lapply(1:indiv_chunk, function(x){round(((df_mom[x,] + df_dad[x,])/2)+rnorm(SNP,mean = 0, sd = 0.005 ))}) %>%
        do.call(rbind, .)
      df_sib <- Normalize_data(df_sib, MAF)
      l_g_sib <- df_sib %*% beta
      l_e_sib <- rnorm(indiv_chunk, 0, sqrt(1-h))
      l_sib <- l_g_sib + l_e_sib
      sib_matrix[,i] <- ifelse(l_sib >= t, 1, 0)
    }



    #NORMALIZING MATRICES
    df_mom <- Normalize_data(df_mom, MAF) #We can overwrite because df_mom no longer needed
    df_dad <- Normalize_data(df_dad, MAF) #We can overwrite because df_dad no longer needed
    df_child_normal <- Normalize_data(df_child, MAF)


    # GENETISK LIABILITY (child)
    l_g_child <- df_child_normal %*% beta

    #Sampling the environmental liability (child)
    l_e_child <- rnorm(indiv_chunk, 0, sqrt(1-h))

    # Finding the full-liability (child)
    l_child <- l_g_child + l_e_child
    # Using the full-liability to determine binary phenotype (child)
    y_child <- ifelse(l_child >= t, 1, 0)

    ######DAD#########
    # GENETISK LIABILITY (dad) (df_dad is already normalized (overwritten))
    l_g_dad <- df_dad %*% beta
    l_e_dad <- rnorm(indiv_chunk, 0, sqrt(1-h))
    l_dad <- l_g_dad + l_e_dad
    y_dad <- ifelse(l_dad >= t, 1, 0)

    ###### MOM ################

    # GENETISK LIABILITY (mom) (df_mom is already normalized (overwritten))
    l_g_mom <- df_mom %*% beta
    l_e_mom <- rnorm(indiv_chunk, 0, sqrt(1-h))
    l_mom <- l_g_mom + l_e_mom
    y_mom <- ifelse(l_mom >= t, 1, 0)



    true_vals <- cbind(id_fra:id_til,l_g_child,l_e_child, l_child, y_child, l_g_dad, l_e_dad, l_dad, y_dad, l_g_mom, l_e_mom, l_mom, y_mom, sib_matrix)
    colnames(true_vals) <- c("ID", "LG_child", "LE_child", "L_child","Child","LG_dad", "LE_dad", "L_dad", "Dad", "LG_mom", "LE_mom", "L_mom", "Mom", "Sib1", "Sib2", "Sib3", "Sib4", "Sib5", "Sib6", "Sib7")



    rm(df_dad,df_mom,df_child_normal, sib_matrix, df_sib) #no longer needed (free space)

    for_plink <- TransformToPlink(df_child,indiv_chunk)
    for_plink <- cbind(id_fra:id_til,id_fra:id_til,for_plink)

    locked = flock::lock(lock) # locks file
    fwrite(as.data.table(true_vals),
           paste("./TRUE","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep=""),
           quote = F,
           sep = " ",
           col.names = F,
           append = T)

    #Filling data in the for_plink matrix (for the first #indiv_chunk persons)
    fwrite(as.data.table(for_plink),
           paste("./DATA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".ped", sep=""),
           quote = F,
           sep = " ",
           append = T)

    # creating and appending to Phone-file
    PHENO <- cbind(id_fra:id_til,id_fra:id_til,y_child)
    fwrite(as.data.table(PHENO),
           paste("./Pheno","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,".txt", sep=""),
           quote = F,
           sep = " ",
           append = T)

    unlock(locked) #unlocks file
    rm(for_plink,PHENO)
  }, future.seed = T)

  writeLines(paste("DONE!"))

}
