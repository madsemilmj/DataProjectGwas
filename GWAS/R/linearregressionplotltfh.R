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
  assoc_file<- paste("./data/LTFH","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h,"_",5,"_",sib,".qassoc", sep="")
  if (file.exists(assoc_file)) {

    assoc_fileltfh <- data.table::fread(assoc_file)
    beta_file <- data.table::fread(paste("./data/BETA","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h,"_",5,".txt", sep=""))
    mf <- data.table::fread(paste("./data/MAF","_",format(total_indiv,scientific = F),"_",format(SNP,scientific = F),"_",h,"_",5,".txt", sep=""))

    est <- assoc_fileltfh %>%
      dplyr::filter(SNP == SNPno)

    true <- beta_file %>%
      dplyr::filter(V1 == SNPno)


    true_beta <- (true$V2- mean(beta_file$V2))/sd(beta_file$V2)
    est$BETA_norm <- (est$BETA - mean(assoc_fileltfh$BETA))/sd(assoc_fileltfh$BETA)
    x <- c(0,1,2)
    x_normal<-(x-(2*as.numeric(mf[1, SNPno, with=FALSE])))/(sqrt(2*as.numeric(mf[1,SNPno, with=FALSE])*(1-as.numeric(mf[1,SNPno, with=FALSE]))))
    est_beta = x_normal *est$BETA_norm

    true_beta <- x_normal*true_beta
    df <- tibble::tibble(x = x_normal, est_beta=est_beta, sd= (est$SE- mean(assoc_fileltfh$SE))/sd(assoc_fileltfh$SE), true_beta=true_beta )
    colors <- c("SE-band" = "red", "Est. Beta" = "steelblue", "True Beta" = "black")
    h <- ggplot2::ggplot(df, aes(x=x_normal)) +
      ggplot2::geom_line(ggplot2::aes(y=est_beta, color = "Est. Beta"))+
      ggplot2::geom_line(ggplot2::aes(y=true_beta, color = "True Beta"))+
      ggplot2::geom_line(ggplot2::aes(y=est_beta-sd, color="SE-band"), linetype="twodash")+
      ggplot2::geom_line(ggplot2::aes(y=est_beta+sd, color="SE-band"), linetype="twodash")+
      ggplot2::labs(x = paste("SNP number: ",SNPno,sep=""), y = "Pheno",color="Legend") +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::theme_light()
    return(h)
  } else {
    print("No data exists! - Try another predefined or run our simulation")
  }

}
