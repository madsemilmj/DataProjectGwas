#' manhattanGWAXSNP plot
#'
#' This function plots a manhatten plot and highlight a specific SNP for the GWAX-method
#' @param SNPno The specific SNP we are looking at
#' @param total_indiv Number of indicviduls in the dataset
#' @param SNP Is the number of SNPs in the data
#' @param h Is the heritability usually 0.5
#' @param sib Is a binary for indicating if we look at sibling history
#' @keywords manhattanGWAXSNP
#' @export
#' @importFrom dplyr %>%
#' @examples
#' manhattanGWAXSNP(SNPno = 10, total_indiv = 1000, SNP = 1000, h = 0.5, sib = 1)

manhattanGWAXSNP <- function(SNPno, total_indiv, SNP, h, sib){
  stringer <- paste("./data/GWAX","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_5_",sib,".assoc", sep="")
  if (file.exists(stringer)) {
    assoc_file <- fread(stringer)
    assoc_file$is_annotate <- ifelse(assoc_file$SNP==SNPno,1,0)
    don <- assoc_file %>%
      dplyr::filter(-log10(P)>1 | is_annotate == 1) # to make plot lighter

    plot1 <- ggplot2::ggplot(don, ggplot2::aes(x=SNP, y=-log10(P))) +

      # Show all points
      ggplot2::geom_point(color="darkgrey", alpha=0.8, size=1) +
      ggplot2::geom_point(data=subset(don, is_annotate==1), color="red", size=2.5) +

      #geom_point(data=subset(don, -log10(P)>thr), color="orange", size=2.5) +
      ggrepel::geom_label_repel(data=subset(don, is_annotate==1), ggplot2::aes(label=SNP),
                                fontface = 'bold',
                                box.padding = 10,
                                point.padding = 0.75,
                                nudge_x = .15,
                                nudge_y = .5,
                                segment.linetype = 1,
                                arrow = arrow(length = unit(0.015, "npc"))
      ) +

      # Custom the theme:
      ggplot2::theme_light() +

      ggplot2::theme(
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
    return(plot1)
  } else {
    print("No data exists! - Try another predefined or run our simulation")
  }
}
