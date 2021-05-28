#' Manhattan Big Plot LT-FH
#'
#' This function plots a manhatten plot with tooltip
#' @param total_indiv Number of indicviduls in the dataset
#' @param SNP Is the number of SNPs in the data
#' @param h Is the heritability usually 0.5
#' @param sib Is a binary for indicating if we look at sibling history
#' @param BFC Is a binary variable indicating if we look at Bon-ferroni-corrected p-values.
#' @keywords Manhattan
#' @export
#' @importFrom dplyr %>%
#' @examples
#' manhattanBigCC(total_indiv = 1000, SNP = 1000, h = 0.5, sib=0, BFC = 1)
manhattanBigltfh <- function(total_indiv, SNP, h, sib=0, BFC = 1){
  file <- paste("./data/LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5_",sib,".qassoc", sep="")
  if (file.exists(file)){
    # Read files
    assoc_string <- paste("./data/LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5_",sib,".qassoc", sep="")
    assoc_file <- data.table::fread(assoc_string)
    causal_string <- paste("./data/BETA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5.txt", sep="")
    causal1 <- data.table::fread(causal_string)
    causal1$C <- ifelse(causal1$V2 != 0, "yes", "no")
    assoc_file <- dplyr::left_join(assoc_file,dplyr::select(causal1,V1,C), by = c("SNP" = "V1"))


    # Threshold
    thr = -log10(.05/BFC)

    don <- assoc_file %>%
      # Filter SNP to make the plot lighter
      dplyr::filter(-log10(P)>1)

    # Prepare text description for each SNP:
    don$text <- paste("SNP: ", don$SNP, "\nLOD score:", -log10(don$P) %>% round(2), "\nCausal SNP:", don$C ,sep="")

    # Make the plot
    p <- ggplot2::ggplot(don, ggplot2::aes(x=SNP, y=-log10(P), text=text)) +
      # Show all points
      ggplot2::geom_point(ggplot2::aes(color=ifelse(C == "yes" & -log10(P)<=thr, 'darkorange1', 'darkgrey')), alpha=0.8, size=1.3) +
      ggplot2::scale_color_manual(values = c("darkorange1" = "darkorange1", "darkgrey" = "darkgrey")) +

      # Add causal points that are also significant
      ggplot2::geom_point(data=subset(don, C != 'no' & -log10(P)>thr), color="cornflowerblue", size=1.6) +

      # Add evt. points that are found significant but not causal (false-positive)
      ggplot2::geom_point(data=subset(don, C == 'no' & -log10(P)>thr), color="red", size=1.6)+

      # Custom the theme:
      ggplot2::theme_light() +

      ggplot2::theme(
        legend.position="none",
        panel.border = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank()
      )
    return(plotly::ggplotly(p, tooltip="text"))
  } else {
    print("No data exists! - Try another predefined or run our simulation")
  }
}
