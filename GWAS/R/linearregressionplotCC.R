#' LinearRegressionPlot (Beta)
#'
#' This function plots linearregresion plot to illustrate beta for case_ctrl
#' @param SNPno The specific SNP we are looking at
#' @param total_indiv Number of indicviduls in the dataset
#' @param SNP Is the number of SNPs in the data
#' @param h Is the heritability usually 0.5
#' @param sib Is a binary for indicating if we look at sibling history
#' @keywords regplot
#' @export
#' @importFrom dplyr %>%
#' @examples
#' linearregressionplotcc(SNPno = 10, total_indiv = 1000, SNP = 1000, h = 0.5, sib=0)
linearregressionplotcc <- function(SNPno, total_indiv, SNP, h, sib){
  assoc_file<- paste("./data/case_ctrl","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h,"_",5,".qassoc", sep="")
  if (file.exists(assoc_file)) {
    assoc_fileCC <- data.table::fread(assoc_file)
    beta_file <- data.table::fread(paste("./data/BETA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h,"_",5,".txt", sep=""))
    mf <- data.table::fread(paste("./data/MAF","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h,"_",5,".txt", sep=""))

    est <- assoc_fileCC %>%
      dplyr::filter(SNP == SNPno)

    true <- beta_file %>%
      dplyr::filter(V1 == SNPno)


    true_beta <- (true$V2- mean(beta_file$V2))/sd(beta_file$V2)
    est$BETA_norm <- (est$BETA - mean(assoc_fileCC$BETA))/sd(assoc_fileCC$BETA)
    x <- c(0,1,2)
    x_normal<-(x-(2*as.numeric(mf[1, SNPno, with=FALSE])))/(sqrt(2*as.numeric(mf[1,SNPno, with=FALSE])*(1-as.numeric(mf[1,SNPno, with=FALSE]))))
    est_beta = x_normal * est$BETA_norm

    true_beta <- x_normal*true_beta

    upperband <- (est$BETA_norm + (est$SE- mean(assoc_fileCC$SE, na.rm = T))/sd(assoc_fileCC$SE, na.rm = T))*x_normal
    lowerband <- (est$BETA_norm - (est$SE- mean(assoc_fileCC$SE, na.rm = T))/sd(assoc_fileCC$SE, na.rm = T))*x_normal

    df <- tibble::tibble(x = x_normal, est_beta=est_beta, upperband= upperband, lowerband = lowerband, true_beta=true_beta )
    colors <- c("SE-band" = "red", "Est. Beta" = "steelblue", "True Beta" = "black")
    h <- ggplot2::ggplot(df, ggplot2::aes(x=x_normal)) +
      ggplot2::geom_line(ggplot2::aes(y=est_beta, color = "Est. Beta"))+
      ggplot2::geom_line(ggplot2::aes(y=true_beta, color = "True Beta"))+
      ggplot2::geom_line(ggplot2::aes(y=lowerband, color="SE-band"), linetype="twodash")+
      ggplot2::geom_line(ggplot2::aes(y=upperband, color="SE-band"), linetype="twodash")+
      ggplot2::labs(x = paste("SNP number: ",SNPno,sep=""), y = "Pheno",color="Legend") +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::theme_light()
    return(h)
  } else {
    print("No data exists! - Try another predefined or run our simulation")
  }

}
