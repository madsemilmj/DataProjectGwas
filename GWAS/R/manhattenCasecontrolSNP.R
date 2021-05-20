#' manhattancasectrlSNP plot
#'
#' This function plots a manhatten plot and highlight a specific SNP for the case-control-method
#' @param SNPno The specific SNP we are looking at
#' @param total_indiv Number of indicviduls in the dataset
#' @param SNP Is the number of SNPs in the data
#' @param h Is the heritability usually 0.5
#' @param sib Is a binary for indicating if we look at sibling history
#' @keywords manhattancasectrlSNP
#' @export
#' @importFrom dplyr %>%
#' @examples
#' manhattancasectrlSNP(SNPno = 10, total_indiv = 1000, SNP = 1000, h = 0.5, sib=0)

manhattancasectrlSNP <- function(SNPno, total_indiv, SNP, h, sib){
  stringer <- paste("./data/case_ctrl","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5.assoc", sep="")
  if (file.exists(stringer)) {
    assoc_file <- data.table::fread(stringer)
    assoc_file$is_annotate <- replicate(nrow(assoc_file),0)
    assoc_file$is_annotate[SNPno] <- 1
    don <- assoc_file %>% 
      dplyr::filter(-log10(P)>1 | SNP == SNPno) # to make plot lighter
    
    plot3 <- ggplot2::ggplot(don, ggplot2::aes(x=SNP, y=-log10(P))) +
      
      # Show all points
      ggplot2::geom_point(color="darkgrey", alpha=0.8, size=1) +
      ggplot2::geom_point(data=subset(don, SNP==SNPno), color="cornflower blue", size=2.5) +
      
      #geom_point(data=subset(don, -log10(P)>thr), color="orange", size=2.5) +
      ggrepel::geom_label_repel(data=subset(don, is_annotate==1), ggplot2::aes(label=SNP), 
                                size=4, box.padding = 0.1, point.padding =1,
                                arrow = arrow(length = unit(0.02, "npc")),ylim=5.5) +
      
      # Custom the theme:
      ggplot2::theme_light() +
      
      ggplot2::theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
    return(plot3)
    
  } else {
    print("No data exists! - Try another predefined or run our simulation")
  }
}