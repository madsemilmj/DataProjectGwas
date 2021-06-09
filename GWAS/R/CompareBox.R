#' Boxplot for comparisons
#'
#' This function plots 3 boxplots to compare the 3 methods
#' @param sib Binary indicator to indicate if we want to look at siblings
#' @keywords Boxplots
#' @export
#' @importFrom dplyr %>%
#' @examples
#' CompareBox(sib = 1, includeTG = 0)
CompareBox <- function(sib, includeTG){
  total_indiv <- 100000
  SNP <- 100000
  h <- 0.5
  file <- paste("./data/BoxData","_",sib,".txt", sep="")
  if (file.exists(file)){
    df <- data.table::fread(file)
    df <- df %>%
      dplyr::filter(Method == "LT-FH" & MeanCausal<38 | Method == "CaseControl" & MeanCausal<38 | Method == "GWAX" & MeanCausal<38 | Method == "TrueGeneticLiab" )
    if (includeTG == 0){
      df <- df %>%
        dplyr::filter(Method != "TrueGeneticLiab")
    }
    NullPlot <- ggplot2::ggplot(df, ggplot2::aes(x=Method, y=MeanNull))+
      ggplot2::geom_boxplot()+
      ggplot2::labs(title = latex2exp::TeX('Average null $\\chi^2$'))+
      ggplot2::theme_light()+
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5))

    CausalPlot <- ggplot2::ggplot(df, ggplot2::aes(x=Method, y=MeanCausal))+
      ggplot2::geom_boxplot()+
      ggplot2::labs(title = latex2exp::TeX('Average causal $\\chi^2$'))+
      ggplot2::theme_light()+
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5))

    PowerPlot <- ggplot2::ggplot(df, ggplot2::aes(x=Method, y=Power))+
      ggplot2::geom_boxplot()+
      ggplot2::labs(title = "Power")+
      ggplot2::theme_light()+
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5))
    return(cowplot::plot_grid(NullPlot,CausalPlot,PowerPlot, ncol = 3))
  } else {
    print("No data exists! - Try another predefined or run our simulation")
  }

}


