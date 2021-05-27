#' ConfusionMatrixLtfh
#'
#' This function plots confusion matrix for the LT-FH method
#' @param total_indiv Number of indicviduls in the dataset
#' @param SNP Is the number of SNPs in the data
#' @param h Is the heritability usually 0.5
#' @param sib Is a binary for indicating if we look at siblings
#' @param BFC Is a variable defining the size of the devider (BFC)
#' @keywords ConfusionMatrix
#' @export
#' @importFrom dplyr %>%
#' @examples
#' ConfusionMatrixLtfh(total_indiv = 1000, SNP = 1000, h = 0.5, sib =1, BFC = 1)
ConfusionMatrixGwax <- function(total_indiv, SNP, h, sib, BFC){
  thr = -log10(.05/BFC)
  #Fetching relevant data
  assoc_string <- paste("./data/LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5_",sib,".qassoc", sep="")
  if (file.exists(assoc_string)){
    assoc_file <- data.table::fread(assoc_string)
    causal_string <- paste("./data/BETA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5.txt", sep="")
    causal1 <- data.table::fread(causal_string)
    #Small data-manipulation
    causal1$C <- ifelse(causal1$V2 != 0, "yes", "no")
    assoc_file <- dplyr::left_join(assoc_file,dplyr::select(causal1,V1,C), by = c("SNP" = "V1"))
    assoc_file$predC <- ifelse(-log10(assoc_file$P)>thr,1,0)

    # Preparing data for confusion matrix
    truth_count_causal <- sum(ifelse(assoc_file$C == "yes",1,0))
    truth_count_Ncausal <- nrow(assoc_file)-truth_count_causal

    pred_count_causal <- sum(assoc_file$predC)
    pred_count_Ncausal <- nrow(assoc_file)-pred_count_causal

    #first col
    caus_caus <- sum(ifelse(assoc_file$C == "yes" & assoc_file$predC == 1,1,0))
    caus_notcaus <- sum(ifelse(assoc_file$C == "yes" & assoc_file$predC == 0,1,0))

    #second_col
    notcaus_caus <- sum(ifelse(assoc_file$C == "no" & assoc_file$predC == 1,1,0))
    notcaus_notcaus <- sum(ifelse(assoc_file$C == "no" & assoc_file$predC == 0,1,0))


    lvs <- c("Causal", "Not Causal")
    truth <- factor(rep(lvs, times = c(truth_count_causal, truth_count_Ncausal)),
                    levels = rev(lvs))

    pred <- factor(
      c(
        rep(lvs, times = c(caus_caus, caus_notcaus)),
        rep(lvs, times = c(notcaus_caus, notcaus_notcaus))),
      levels = rev(lvs))


    #Preparing data for plotting
    table <- data.frame(caret::confusionMatrix(pred, truth, dnn = c("Test", "True"))$table)

    plotTable <- table %>%
      dplyr::mutate(Label = ifelse(table$Test == table$True, "Correct", "Not_correct")) %>%
      dplyr::group_by(True)

    #Plot
    p <- ggplot2::ggplot(data = plotTable, mapping = ggplot2::aes(x = True, y = Test, fill = Label)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1, size = 14) +
      ggplot2::scale_fill_manual(values = c(Correct = "green", Not_correct = "red")) +
      ggplot2::theme_light() +
      ggplot2::xlim(rev(levels(table$True)))
    return(p)


  } else {
    print("No data exists for this ConfusionMatrixLtfh! - Try another predefined or run our simulation")
  }
}
