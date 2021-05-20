#' LinearRegressionPlot (Beta)
#'
#' This function plots linearregresion plot to illustrate beta
#' @param SNPno The specific SNP we are looking at
#' @param total_indiv Number of indicviduls in the dataset
#' @param SNP Is the number of SNPs in the data
#' @param h Is the heritability usually 0.5
#' @param sib Is a binary for indicating if we look at sibling history
#' @keywords DistPlot
#' @export
#' @importFrom dplyr %>%
#' @examples
#' linearregressionplotltfh(SNPno = 10, total_indiv = 1000, SNP = 1000, h = 0.5, sib=0)


linearregressionplotltfh <- function(SNPno, total_indiv, SNP, h, sib){
  assoc_file<- paste("./data/LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h*100,"_",5,"_",sib,".qassoc", sep="")
  if (file.exists(assoc_file)) {
    assoc_fileltfh <- data.table::fread(assoc_file)
    beta_file <- data.table::fread(paste("./data/BETA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h,"_",5,".txt", sep=""))
    #Creating bonferroni adjusted p-values
    est <- assoc_fileltfh %>%
      dplyr::filter(SNP == SNPno)

    true <- beta_file %>%
      dplyr::filter(V1 == SNPno)

    true_beta <- true$V2

    df <- tibble::tibble(x = c(0,1), est_beta = 0:1 * est$BETA, sd = est$SE, true_beta = 0:1*true_beta)
    colors <- c("SE-band" = "red", "Est. Beta" = "steelblue", "True Beta" = "black")
    h <- ggplot2::ggplot(df, ggplot2::aes(x=x)) +
      ggplot2::geom_line(aes(y=est_beta, color = "Est. Beta"))+
      ggplot2::geom_line(aes(y=true_beta, color = "True Beta"))+
      ggplot2::geom_line(aes(y=est_beta-sd, color="SE-band"), linetype="twodash")+
      ggplot2::geom_line(aes(y=est_beta+sd, color="SE-band"), linetype="twodash")+
      ggplot2::labs(x = paste("SNP number: ",SNPno,sep=""), y = "Pheno",color="Legend") +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::theme_light()
    return(h)
  } else {
      print("No data exists! - Try another predefined or run our simulation")
    }

}
