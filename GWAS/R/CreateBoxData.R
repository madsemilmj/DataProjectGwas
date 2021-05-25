#' Creating data for boxplot
#'
#' This function creates data used to plot 3 boxplots to compare the 3 methods
#' @param sib Binary indicator to indicate if we want to look at siblings
#' @keywords BoxplotsData
#' @export
#' @importFrom dplyr %>%
#' @examples
#' CreateBoxData(sib = 1)
CreateBoxData <- function(sib){
  total_indiv <- 100000
  SNP <- 100000
  h <- 0.5
  file <- paste("./data/LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5_",sib,".qassoc", sep="")
  if (file.exists(file)){
    # Stringers for the first assoc-files
    first_ltfh <- paste("./data/LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5_",sib,".qassoc", sep="")
    first_cc <- paste("./data/case_ctrl","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5.assoc", sep="")
    first_gwax <- paste("./data/GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5_",sib,".assoc", sep="")
    first_TG <- paste("./data/TG","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5.qassoc", sep="")
    #Reading assoc files
    ltfh <- data.table::fread(first_ltfh)
    gwax <- data.table::fread(first_gwax)
    cc <- data.table::fread(first_cc)
    TG <- data.table::fread(first_TG)
    #adding chisq to ltfhand TG
    ltfh$CHISQ <- qchisq(ltfh$P,df=1,lower.tail = FALSE)
    TG$CHISQ <- qchisq(TG$P,df=1,lower.tail = FALSE)
    ##SLET
    cc$CHISQ <- qchisq(cc$P,df=1,lower.tail = FALSE)
    gwax$CHISQ <- qchisq(gwax$P,df=1,lower.tail = FALSE)
    ##SLET
    #Check for found causals
    ltfh$causal <- ifelse(ltfh$P < (0.05)/1000000,1,0)
    gwax$causal <- ifelse(gwax$P < (0.05)/1000000,1,0)
    cc$causal <- ifelse(cc$P < (0.05)/1000000,1,0)
    TG$causal <- ifelse(TG$P < (0.05)/1000000,1,0)
    #calculating-power
    ltfhPower <- nrow(subset(ltfh, causal == 1))/1000
    gwaxPower <- nrow(subset(gwax, causal == 1))/1000
    ccPower <- nrow(subset(cc, causal == 1))/1000
    TGPower <- nrow(subset(TG, causal == 1))/1000
    #stringer for true-beta
    beta_string <-  paste("./data/BETA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5.txt", sep="")
    true_beta <- data.table::fread(beta_string)
    # Joining to get true causals
    causals <- true_beta %>%
      dplyr::filter(V2 != 0) %>%
      dplyr::left_join(dplyr::select(ltfh,SNP,CHISQ), by = c("V1" = "SNP")) %>%
      dplyr::rename("SNP" = "V1", "TrueBeta" = "V2", "CHISQ_ltfh" = "CHISQ") %>%
      dplyr::left_join(dplyr::select(gwax,SNP,CHISQ), by = c("SNP" = "SNP")) %>%
      dplyr::rename("CHISQ_gwax" = "CHISQ") %>%
      dplyr::left_join(dplyr::select(cc,SNP,CHISQ), by = c("SNP" = "SNP")) %>%
      dplyr::rename("CHISQ_cc" = "CHISQ") %>%
      dplyr::left_join(dplyr::select(TG,SNP,CHISQ), by = c("SNP" = "SNP")) %>%
      dplyr::rename("CHISQ_TG" = "CHISQ")

    #joining to get nulls
    nulls <- true_beta %>%
      dplyr::filter(V2 == 0) %>%
      dplyr::left_join(dplyr::select(ltfh,SNP,CHISQ), by = c("V1" = "SNP")) %>%
      dplyr::rename("SNP" = "V1", "TrueBeta" = "V2", "CHISQ_ltfh" = "CHISQ") %>%
      dplyr::left_join(dplyr::select(gwax,SNP,CHISQ), by = c("SNP" = "SNP")) %>%
      dplyr::rename("CHISQ_gwax" = "CHISQ") %>%
      dplyr::left_join(dplyr::select(cc,SNP,CHISQ), by = c("SNP" = "SNP")) %>%
      dplyr::rename("CHISQ_cc" = "CHISQ") %>%
      dplyr::left_join(dplyr::select(TG,SNP,CHISQ), by = c("SNP" = "SNP")) %>%
      dplyr::rename("CHISQ_TG" = "CHISQ")


    df <- tibble::tibble(MeanNull = mean(nulls$CHISQ_ltfh), MeanCausal = mean(causals$CHISQ_ltfh), Power = ltfhPower, Method = "LT-FH") %>%
      dplyr::add_row(MeanNull = mean(nulls$CHISQ_gwax), MeanCausal = mean(causals$CHISQ_gwax), Power = gwaxPower, Method = "GWAX") %>%
      dplyr::add_row(MeanNull = mean(nulls$CHISQ_cc), MeanCausal = mean(causals$CHISQ_cc), Power = ccPower, Method = "CaseControl") %>%
      dplyr::add_row(MeanNull = mean(nulls$CHISQ_TG), MeanCausal = mean(causals$CHISQ_TG), Power = TGPower, Method = "TrueGeneticLiab")

    #Looping through rest of files (if they exsist) and add them
    for (i in 1:9){
      file <- paste("./data/NineSims/",i,"/data/LTFH_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5_",sib,".qassoc", sep="")
      if (file.exists(file)){
        # Stringers for the current assoc-files
        ltfh_string <- paste("./data/NineSims/",i,"/data/LTFH_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5_",sib,".qassoc", sep="")
        cc_string <- paste("./data/NineSims/",i,"/data/case_ctrl","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5.assoc", sep="")
        gwax_string <- paste("./data/NineSims/",i,"/data/GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5_",sib,".assoc", sep="")
        TG_string <- paste("./data/NineSims/",i,"/data/TG","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5.qassoc", sep="")
        #Reading assoc files
        ltfh <- data.table::fread(ltfh_string)
        gwax <- data.table::fread(gwax_string)
        cc <- data.table::fread(cc_string)
        TG <- data.table::fread(TG_string)
        #adding chisq to ltfh and TG
        ltfh$CHISQ <- qchisq(ltfh$P,df=1,lower.tail = FALSE)
        TG$CHISQ <- qchisq(TG$P,df=1,lower.tail = FALSE)
        ##SLET
        cc$CHISQ <- qchisq(cc$P,df=1,lower.tail = FALSE)
        gwax$CHISQ <- qchisq(gwax$P,df=1,lower.tail = FALSE)
        ##SLET
        #Check for found causals
        ltfh$causal <- ifelse(ltfh$P < (0.05)/1000000,1,0)
        gwax$causal <- ifelse(gwax$P < (0.05)/1000000,1,0)
        cc$causal <- ifelse(cc$P < (0.05)/1000000,1,0)
        TG$causal <- ifelse(TG$P < (0.05)/1000000,1,0)
        #calculating-power
        ltfhPower <- nrow(subset(ltfh, causal == 1))/1000
        gwaxPower <- nrow(subset(gwax, causal == 1))/1000
        ccPower <- nrow(subset(cc, causal == 1))/1000
        TGPower <- nrow(subset(TG, causal == 1))/1000
        #stringer for true-beta
        beta_string <-  paste("./data/NineSims/",i,"/data/BETA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5.txt", sep="")
        true_beta <- data.table::fread(beta_string)
        # Joining to get true causals
        causals <- true_beta %>%
          dplyr::filter(V2 != 0) %>%
          dplyr::left_join(dplyr::select(ltfh,SNP,CHISQ), by = c("V1" = "SNP")) %>%
          dplyr::rename("SNP" = "V1", "TrueBeta" = "V2", "CHISQ_ltfh" = "CHISQ") %>%
          dplyr::left_join(dplyr::select(gwax,SNP,CHISQ), by = c("SNP" = "SNP")) %>%
          dplyr::rename("CHISQ_gwax" = "CHISQ") %>%
          dplyr::left_join(dplyr::select(cc,SNP,CHISQ), by = c("SNP" = "SNP")) %>%
          dplyr::rename("CHISQ_cc" = "CHISQ") %>%
          dplyr::left_join(dplyr::select(TG,SNP,CHISQ), by = c("SNP" = "SNP")) %>%
          dplyr::rename("CHISQ_TG" = "CHISQ")


        #joining to get nulls
        nulls <- true_beta %>%
          dplyr::filter(V2 == 0) %>%
          dplyr::left_join(dplyr::select(ltfh,SNP,CHISQ), by = c("V1" = "SNP")) %>%
          dplyr::rename("SNP" = "V1", "TrueBeta" = "V2", "CHISQ_ltfh" = "CHISQ") %>%
          dplyr::left_join(dplyr::select(gwax,SNP,CHISQ), by = c("SNP" = "SNP")) %>%
          dplyr::rename("CHISQ_gwax" = "CHISQ") %>%
          dplyr::left_join(dplyr::select(cc,SNP,CHISQ), by = c("SNP" = "SNP")) %>%
          dplyr::rename("CHISQ_cc" = "CHISQ") %>%
          dplyr::left_join(dplyr::select(TG,SNP,CHISQ), by = c("SNP" = "SNP")) %>%
          dplyr::rename("CHISQ_TG" = "CHISQ")

        df <- df %>%
          dplyr::add_row(MeanNull = mean(nulls$CHISQ_ltfh), MeanCausal = mean(causals$CHISQ_ltfh), Power = ltfhPower, Method = "LT-FH") %>%
          dplyr::add_row(MeanNull = mean(nulls$CHISQ_gwax), MeanCausal = mean(causals$CHISQ_gwax), Power = gwaxPower, Method = "GWAX") %>%
          dplyr::add_row(MeanNull = mean(nulls$CHISQ_cc), MeanCausal = mean(causals$CHISQ_cc), Power = ccPower, Method = "CaseControl") %>%
          dplyr::add_row(MeanNull = mean(nulls$CHISQ_TG), MeanCausal = mean(causals$CHISQ_TG), Power = TGPower, Method = "TrueGeneticLiab")


      }else {print("no file")}
    }

    data.table::fwrite(data.table::as.data.table(df),
           paste("./data/BoxData","_",sib,".txt", sep=""),
           quote = F,
           sep = " ",
           col.names = T)

  } else {
    print("No data exists! - Try another predefined or run our simulation")
  }

}
