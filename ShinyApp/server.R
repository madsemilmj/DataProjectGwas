#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(devtools)

library(shiny)
library(knitr)
library(markdown)
library(DT)

library(tidyverse)
library(knitr)
library(data.table)
library(vroom)
library(future) 
library(flock)
library(cowplot)
library(ggplot2)
library(plotly)
library(dplyr)
library(shinythemes)
library(GWAS)
library(purrr)
library(stats)
library(graphics)  
library(mvtnorm)



get_assoc_cc <- function(total_indiv, SNP, h, k){
  assoc_file<- paste("./data/case_ctrl","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h,"_",k*100,".qassoc", sep="")
  (fread(assoc_file))
  
}

get_assoc_ltfh <- function(total_indiv, SNP, h, k, sib){
  assoc_file<- paste("./data/LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h,"_",k*100,"_",sib,".qassoc", sep="")
  (fread(assoc_file))}

get_assoc_GW <- function(total_indiv, SNP, h, k, sib){
  assoc_file<- paste("./data/GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h,"_",k*100,"_",sib,".qassoc", sep="")
  (fread(assoc_file))
}


choose1 <- function(x){
  if (length(x)==0){
  return(0)
  }else if (grepl(x[1], "Child is a case")){
    return(1)
  }else {
    return(0)
  }}


choose2 <- function(x){
  if (length(x)==0){
    return(0)
  }else if (grepl(x[1], "Mother is a case")){
    return(1)
  }else if (length(x)==1){
      return(0)
  }else if (grepl(x[2], "Mother is a case")){
    return(1)
    }else {
    return(0)
  }}

choose3 <- function(x){
  if (length(x)==0){
    return(0)
  }else if (grepl(x[1], "Father is a case")){
    return(1)
  }else if (length(x)==1){
    return(0)
  }else if (grepl(x[2], "Father is a case")){
    return(1)
  }else if (length(x)==2){
    return(0)
  }else if (grepl(x[3], "Father is a case")){
    return(1)
    }else {
    return(0)
  }}











# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    output$Predifined <- renderPrint({ input$predifined })
    
    
    
    output$upload <- renderPrint( {str(input$file)})
    
  
    output$upload1 <- renderPrint({
      str(input$file1)
    })
      
    output$upload2 <- renderPrint({
      str(input$file2)
    })
    
    output$upload3 <- renderPrint({
      str(input$file3)
    })

    output$uploader <-reactive({if(is.null(input$upload)){
      "Please upload some data"
    }else{
      "The data format is not correct"
    }})
    
    output$uploader1  <-reactive({if(is.null(input$upload1)){
      "Please upload some data"
    }else{
      "The data format is not correct"
    }})
    output$uploader2  <-reactive({if(is.null(input$upload2)){
      "Please upload some data"
    }else{
      "The data format is not correct"
    }})
    output$uploader3  <-reactive({if(is.null(input$upload3)){
      "Please upload some data"
    }else{
      "The data format is not correct"
    }})
  
    

    
    simulated <- observeEvent(input$run_sim, {MasterFunc(total_indiv = input$Individuals_sim, indiv_chunk = 100, SNP =  input$SNPs_sim, h= input$H2_sim, c=input$SNPs_sim/100, k = 0.05)
    output$simulator <- renderPrint({"Done!"}) 
   })
      
  

    output$casectrlsim <- renderDT({get_assoc_cc(input$Individuals_sim, input$SNPs_sim,input$H2_sim*100, k = 0.05 ) }, selection = 'single')
    
    output$GWAXheadsim <- renderDT({get_assoc_GW(input$Individuals_sim, input$SNPs_sim,input$H2_sim*100, 0.05, input$Sib_GW_sim)}, selection = 'single') #change sib
    
    
    output$LTFHheadsim <- renderDT({get_assoc_ltfh(input$Individuals_sim, input$SNPs_sim,input$H2_sim*100, 0.05, input$Sib_ltfh_sim )}, selection = 'single') #change sib
    
    simValues <- reactive({
      
      data.frame(
        Name = c("Individuals", "SNPS", "Heritability"),
        Value = as.character(c(input$Individuals_sim, input$SNPs_sim,input$H2_sim)),
        stringsAsFactors = FALSE)
   
    }) 
     output$table_of_sim <- renderTable({
      simValues()}
    )
       
     
     output$table_of_sim1 <- renderTable({
       simValues()}
       
     )
     
     output$table_of_sim2 <- renderTable({
       simValues()}
       
     )
     
     output$table_of_sim3 <- renderTable({
       simValues()}
       
     )
     
     output$table_of_sim4 <- renderTable({
       simValues()}
       
     )

  
    output$casectrl <-renderDT({get_assoc_cc(input$Individuals_cc, input$SNPs_cc,input$h2_cc, 0.05 )}, selection = 'single')
    
    
   
    output$GWAXhead <- renderDT({get_assoc_GW(input$Individuals_GW, input$SNPs_GW,input$h2_GW, 0.05, input$Sib_GW )}, selection = 'single')
     
    output$LTFHhead <- renderDT({get_assoc_ltfh(input$Individuals_lfth, input$SNPs_ltfh,input$h2_ltfh, 0.05, input$Sib_ltfh )}, selection = 'single')
    


    
    output$plotsteriot_pre <- renderPlot({plosteriot(as.numeric(input$Individuals_lfth), as.numeric(input$SNPs_ltfh) ,as.numeric(input$h2_ltfh)/100, choose1(input$Plotsterior), choose2(input$Plotsterior), choose3(input$Plotsterior),  input$Sibs,input$Sibs_case )})
    
    
    output$plotsteriot_sim <- renderPlot({plosteriot(as.numeric(input$Individuals_sim), as.numeric(input$SNPs_sim) ,as.numeric(input$H2_sim)/100, choose1(input$Plotsterior), choose2(input$Plotsterior), choose3(input$Plotsterior),  input$Sibs,input$Sibs_case )    })
    
    

    
    output$RegressionCCsim<- renderPlot({
      if (length(input$casectrlsim_rows_selected)){
        linearregressionplotcc(ifelse(input$casectrlsim_rows_selected==1, input$casectrlsim_rows_selected, (input$casectrlsim_rows_selected-1)) ,as.numeric(input$Individuals_sim), as.numeric(input$SNPs_sim), as.numeric(input$H2_sim)*100, as.numeric(input$sib_sim )) }
    }) 
    
    
    output$manhattanCCsim <- renderPlot({
      if (length(input$casectrlsim_rows_selected)){
        manhattancasectrlSNP(ifelse(input$casectrlsim_rows_selected==1, input$casectrlsim_rows_selected, (input$casectrlsim_rows_selected-1)),input$Individuals_sim, input$SNPs_sim, as.numeric(input$H2_sim), input$sib_sim ) }
    })
    
    
    output$RegressionCC <- renderPlot({
      if (length(input$casectrl_rows_selected)){
        linearregressionplotcc(ifelse(input$casectrl_rows_selected==1, input$casectrl_rows_selected, (input$casectrl_rows_selected-1)),input$Individuals_cc, input$SNPs_cc, as.numeric(input$h2_cc), input$Sib_cc ) }
    }) 
    
    
    output$manhattanCC <- renderPlot({
      if (length(input$casectrl_rows_selected)){
        manhattancasectrlSNP(ifelse(input$casectrl_rows_selected==1, input$casectrl_rows_selected, (input$casectrl_rows_selected-1)),input$Individuals_cc, input$SNPs_cc, as.numeric(input$h2_cc)/100, input$Sib_cc ) }
    })
    
    
    
    output$RegressionGWAXsim<- renderPlot({
      if (length(input$GWAXheadsim_rows_selected)){
        linearregressionplotgwax(ifelse(input$GWAXheadsim_rows_selected==1, input$GWAXheadsim_rows_selected, (input$GWAXheadsim_rows_selected-1)),as.numeric(input$Individuals_sim), as.numeric(input$SNPs_sim), as.numeric(input$H2_sim)*100, as.numeric(input$sib_sim )) }
    }) 
    
    
    output$manhattanGWAXsim <- renderPlot({
      if (length(input$GWAXheadsim_rows_selected)){
        manhattanGWAXSNP(ifelse(input$GWAXheadsim_rows_selected==1, input$GWAXheadsim_rows_selected, (input$GWAXheadsim_rows_selected-1)),input$Individuals_sim, input$SNPs_sim, as.numeric(input$H2_sim), input$sib_sim ) }
    })
    
    
    output$RegressionGWAX <- renderPlot({
      if (length(input$GWAXhead_rows_selected)){
        linearregressionplotgwax(ifelse(input$GWAXhead_rows_selected==1, input$GWAXhead_rows_selected, (input$GWAXhead_rows_selected-1)),input$Individuals_GW, input$SNPs_GW,input$h2_GW, input$Sib_GW ) }
    }) 
    
    
    output$manhattanGWAX <- renderPlot({
      if (length(input$GWAXhead_rows_selected)){
        manhattanGWAXSNP(ifelse(input$GWAXhead_rows_selected==1, input$GWAXhead_rows_selected, (input$GWAXhead_rows_selected-1)),input$Individuals_GW, input$SNPs_GW, as.numeric(input$h2_GW)/100, input$Sib_GW ) }
    }) 
    
    
    
    
    
    output$Regressionltfhsim<- renderPlot({
      if (length(input$LTFHheadsim_rows_selected)){
        linearregressionplotltfh(ifelse(input$LTFHheadsim_rows_selected==1, input$LTFHheadsim_rows_selected, (input$LTFHheadsim_rows_selected-1)),as.numeric(input$Individuals_sim), as.numeric(input$SNPs_sim), as.numeric(input$H2_sim)*100, as.numeric(input$sib_sim )) }
    }) 
    
    
    output$manhattanLTFHsim <- renderPlot({
      if (length(input$LTFHheadsim_rows_selected)){
        manhattanLTFHSNP(ifelse(input$LTFHheadsim_rows_selected==1, input$LTFHheadsim_rows_selected, (input$LTFHheadsim_rows_selected-1)),input$Individuals_sim, input$SNPs_sim, as.numeric(input$H2_sim), input$sib_sim ) }
    })
    
    
    output$Regressionltfh <- renderPlot({
      if (length(input$LTFHhead_rows_selected)){
     linearregressionplotltfh(ifelse(input$LTFHhead_rows_selected==1, input$LTFHhead_rows_selected, (input$LTFHhead_rows_selected-1)),input$Individuals_lfth, input$SNPs_ltfh,input$h2_ltfh, input$Sib_ltfh ) }
    }) 
  
    
    output$manhattanLTFH <- renderPlot({
      if (length(input$LTFHhead_rows_selected)){
        manhattanLTFHSNP(ifelse(input$LTFHhead_rows_selected==1, input$LTFHhead_rows_selected, (input$LTFHhead_rows_selected-1)),input$Individuals_lfth, input$SNPs_ltfh,as.numeric(input$h2_ltfh)/100, input$Sib_ltfh ) }
    }) 
    
    output$boxcompare <- renderPlot({CompareBox(input$Sib_box1, input$Truegen)})
    
    output$boxcompare1 <- renderPlot({CompareBox(input$Sib_box1, input$Truegen1)})
    output$zscoreGG <- renderPlot({zscorePlotGG(input$Individuals_co, input$SNPs_co, as.numeric(input$h2_co)/100,input$Sib_box1, input$Truegen)      })
    
    output$zscoreGG_sim <-renderPlot({zscorePlotGG(input$Individuals_sim, input$SNPs_sim, as.numeric(input$H2_sim), input$Sib_z_sim, input$Truegen1 )      })

     
    output$bigmanhattancc <-renderPlotly({
      if (input$Choose_type == 1){
        manhattanBigCC(input$Individuals_ma, input$SNPs_ma, as.numeric(input$h2_ma)/100,as.numeric(input$Bonferronicoefficient))
      } else if (input$Choose_type == 2){
        manhattanBiggwax(input$Individuals_ma, input$SNPs_ma, as.numeric(input$h2_ma)/100,input$Sib_box,as.numeric(input$Bonferronicoefficient))
      } else if (input$Choose_type ==3) {
        manhattanBigltfh(input$Individuals_ma, input$SNPs_ma, as.numeric(input$h2_ma)/100,input$Sib_box,as.numeric(input$Bonferronicoefficient))
      } else{
        manhattanBigTG(input$Individuals_ma, input$SNPs_ma, as.numeric(input$h2_ma)/100,as.numeric(input$Bonferronicoefficient))
      }
    })
    
    
    output$bigmanhattansim <-renderPlotly({
      if (input$Choose_type == 1){
        manhattanBigCC(input$Individuals_sim, input$SNPs_sim, as.numeric(input$H2_sim),as.numeric(input$Bonferronicoefficientsim))
      } else if (input$Choose_type == 2){
        manhattanBiggwax(input$Individuals_sim, input$SNPs_sim, as.numeric(input$H2_sim), input$Sib_ma_sim,as.numeric(input$Bonferronicoefficientsim))
      } else if (input$Choose_type == 3){
        anhattanBigltfh(input$Individuals_sim, input$SNPs_sim, as.numeric(input$H2_sim), input$Sib_ma_sim,as.numeric(input$Bonferronicoefficientsim))
      } else {
        manhattanBigTG(input$Individuals_sim, input$SNPs_sim, as.numeric(input$H2_sim)
                       , as.numeric(input$Bonferronicoefficientsim))
      }
    })
    
    output$ConfussionM <- renderPlot({ if (input$Choose_type == 1){
      ConfusionMatrixCC(input$Individuals_ma, input$SNPs_ma, as.numeric(input$h2_ma)/100,as.numeric(input$Bonferronicoefficient))
    } else if (input$Choose_type == 2){
      ConfusionMatrixGwax(input$Individuals_ma, input$SNPs_ma, as.numeric(input$h2_ma)/100,input$Sib_box,as.numeric(input$Bonferronicoefficient))
    } else if (input$Choose_type ==3) {
      ConfusionMatrixLtfh(input$Individuals_ma, input$SNPs_ma, as.numeric(input$h2_ma)/100,input$Sib_box,as.numeric(input$Bonferronicoefficient))
    } else{
      ConfusionMatrixTG(input$Individuals_ma, input$SNPs_ma, as.numeric(input$h2_ma)/100,as.numeric(input$Bonferronicoefficient))
    }
      
    })
    
    
    
    output$ConfussionMsim <- renderPlot({ if (input$Choose_type == 1){
      ConfusionMatrixCC(input$Individuals_sim, input$SNPs_sim, as.numeric(input$H2_sim),as.numeric(input$Bonferronicoefficientsim))
    } else if (input$Choose_type == 2){
      ConfusionMatrixGwax(input$Individuals_sim, input$SNPs_sim, as.numeric(input$H2_sim), input$Sib_ma_sim,as.numeric(input$Bonferronicoefficientsim))
    } else if (input$Choose_type ==3) {
      ConfusionMatrixLtfh(input$Individuals_sim, input$SNPs_sim, as.numeric(input$H2_sim), input$Sib_ma_sim,as.numeric(input$Bonferronicoefficientsim))
    } else{
      ConfusionMatrixTG(input$Individuals_sim, input$SNPs_sim, as.numeric(input$H2_sim) , as.numeric(input$Bonferronicoefficientsim))
    }
      
    })
    
    output$sessionInfo <- renderPrint({
      capture.output(sessionInfo())
    })
    
    output$gibbscode1 <- renderPrint({cat("set.seed(18)
    

h2 <- 0.5
nr_sib <- 1
valT <- 0.05 #prevalence (used to define bounds/case - threshold)

CovMatrix <- GWAS::LinearTransformation(GWAS::CovarianceMatrix(h2,nr_sib))
#Below is where we change the bounds depeding on the configuration
LowerBound <- rep(-Inf,nrow(CovMatrix)-1)
UpperBound <- rep(Inf,nrow(CovMatrix)-1)

CondInfo <- GWAS::ConditionalDistribution(CovMatrix)
sigma <- CondInfo$sigma
mu <- CondInfo$mu
# Total number of samples
S <- 6000
my_sample <- matrix(ncol = (length(LowerBound)+1), nrow = S)
#Below is the initial values
vals <- rep(10,(length(LowerBound)+1))
#Running the gibbs-sampling
for (i in 1:S){
  #Child gen
  vals[1] <- rnorm(1,mu[1,] %*% vals[2:length(vals)], sqrt(sigma[1]))
  #REST
  for (p in 2:(length(LowerBound)+1)){
    mu_prod <- mu[p,] %*% vals[-p]
    l1 <- pnorm(LowerBound[p-1], mu_prod, sqrt(sigma[p]))
    u1 <- pnorm(UpperBound[p-1], mu_prod, sqrt(sigma[p]))
    g1 <- runif(1,l1,u1)
    vals[p] <- qnorm(g1,mu_prod, sqrt(sigma[p]))
  }
  my_sample[i,] <- vals
}")})
    
    

    
output$KODE1 <-  renderPrint({cat("Normalize_data <- function(df,MAF){
    df_normal <- sweep(sweep(df, 2, 2*MAF, '-'), 2, sqrt(2*MAF*(1-MAF)), '/')
    return(df_normal)  }")})   

output$KODE2 <- renderPrint({cat("TransformToPlink <- function(df,indiv_chunk){
  a <- matrix(ncol=3,nrow=2) #C er den riskogivne dvs CC = 2 og AA = 0
  a[1,]<-c('A', 'C', 'C')
  a[2,]<- c('A', 'A', 'C')

  for_plink <- lapply(1:indiv_chunk, function(x){
    tmp <- a[,df[x,]+1]
    c(rep(1,4),rbind(tmp[1,],tmp[2,]))
    })%>%
      do.call(rbind, .)
  return(for_plink)}")})

output$KODE3 <- renderPrint({cat("SimulerData <- function(total_indiv, indiv_chunk, SNP, h, c, k, nr_workers) {")})    
    
    
    
    
    
output$KODE4 <- renderPrint({cat("t <- qnorm(1-k, 0, 1) ")})    
    
    
output$KODE5 <- renderPrint({cat("MAF <- runif(n = SNP, min = 0.01, max = 0.49) # Simulating #SNP probabilities. ")})
    
    
output$KODE6 <- renderPrint({cat('#Generating MAP file for plink
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
         col.names = F)')})    
    
output$KODE7 <- renderPrint({cat("# Drawing C indices from 1:SNP without replacement (represents the C causal SNPS)
  disease <- sample(SNP, c, replace = FALSE, prob = NULL)

  # Using the above indices to find the effect sizes
  beta <- numeric(SNP)
  beta[disease] <- rnorm(c, 0, sqrt(h/c)) ")})   




output$KODE8 <- renderPrint({cat("# Defining genotypes for mom, dad and child
  df_mom <- t(replicate(indiv_chunk, rbinom(SNP, 2, MAF)))
  df_dad <- t(replicate(indiv_chunk, rbinom(SNP, 2, MAF)))")})
    
    
    
    
    
output$KODE9 <- renderPrint({cat("df_child <- lapply(1:indiv_chunk, function(x){round(((df_mom[x,] + df_dad[x,])/2)+rnorm(SNP,mean = 0, sd = 0.005 ))}) %>%
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
  } ")})    
    
    
output$KODE10 <- renderPrint({cat("#NORMALIZING MATRICES
  df_mom <- Normalize_data(df_mom, MAF) #We can overwrite because df_mom no longer needed
  df_dad <- Normalize_data(df_dad, MAF) #We can overwrite because df_dad no longer needed
  df_child_normal <- Normalize_data(df_child, MAF) ")})    
    
    
output$KODE11 <- renderPrint({cat("### DEFINING LIABILITIES
  l_g_child <- df_child_normal %*% beta ")})    
    


output$KODE12 <- renderPrint({cat("#Sampling the environmental liability (child)
  l_e_child <- rnorm(indiv_chunk, 0, sqrt(1-h)) ")})    
    
    
output$KODE13 <- renderPrint({cat("
                                  # Finding the full-liability (child)
  l_child <- l_g_child + l_e_child
  # Using the full-liability to determine binary phenotype (child)
  y_child <- ifelse(l_child >= t, 1, 0) ")})    
    


output$KODE14 <- renderPrint({cat("######DAD#########
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
  y_mom <- ifelse(l_mom >= t, 1, 0)   ")})    


output$KODE15 <- renderPrint({cat('true_vals <- cbind(1:indiv_chunk,l_g_child,l_e_child, l_child, y_child, l_g_dad, l_e_dad, l_dad, y_dad, l_g_mom, l_e_mom, l_mom, y_mom, sib_matrix)
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
         col.names = T) ')})

output$KODE16 <- renderPrint({cat('rm(for_plink,PHENO)
  lock = tempfile()

  #### you generally do not want to put plan inside a function like this. future_lapply collapses to a normal lapply function if no plan is given.
  #### this allows the end user to use it however they want.
  # Defining the number of cores in the processor to run simultaneous
  x <<- nr_workers #GLOBAL VARIABLE SO IT CAN BE USED BY FUTURE
  plan(future::multisession(workers = x))
  nr_it <- total_indiv %/% indiv_chunk #Integer division
  #Performing the steps in parallel
  future.apply::future_lapply(2:nr_it, function(n) {

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

  #writeLines(paste("DONE!"))
  }



master_func <- function(total_indiv, indiv_chunk, SNP, h, c, k){
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
  fwrite(as.data.table(ltfh_pheno_1),
         paste("./Pheno_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_1",".txt", sep=""),
         quote = F,
         sep = " ",
         col.names = T)
  #Making a pheno_file for LTFH without sib
  ltfh_pheno_0 <- Assign_LTFH(Pheno_data = ny_true,
                            valT = k, h2 = h, with_sib = 0)
  fwrite(as.data.table(ltfh_pheno_0),
         paste("./Pheno_LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_0",".txt", sep=""),
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


start <- Sys.time()
master_func(total_indiv = 1000,indiv_chunk = 100, SNP = 300, h = .5, c = 8, k = 0.05)
end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
total_time')})



output$GWASplink <- renderPrint({
  cat('case_ctrl <- paste("case_ctrl","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100, sep="")
shell(cmd =paste("plink --bfile", bfile, "--pheno", phenotype_GWAS, "--pheno-name Pheno --out", case_ctrl, "--assoc")')
})


output$GWAXplink <- renderPrint({
  cat('GWAX <- paste("GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_sib", sep="")
shell(cmd =paste("plink --bfile", bfile, "--pheno", phenotype_GWAX, "--pheno-name Pheno --out", GWAX, "--assoc")')
})


output$LTFHplink <- renderPrint({
  cat(' LTFH <- paste("LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",k*100,"_sib", sep="")
shell(cmd paste("plink --bfile", bfile, "--pheno", phenotype_LTFH, "--pheno-name Pheno --out", LTFH, "--assoc")')
})  
    
    
}) #server end


