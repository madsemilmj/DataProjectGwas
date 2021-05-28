#' Plot distributions
#'
#' This function plots posterior distribution and prior
#' @param total_indiv Number of indicviduls in the dataset
#' @param SNP Is the number of SNPs in the data
#' @param h Is the heritability usually 0.5
#' @param Child Is a binary for indicating case/notcase for child
#' @param Mom Is a binary for indicating case/notcase for Mom
#' @param Dad Is a binary for indicating case/notcase for Dad
#' @param Nr_sib indiacates the number of siblings for the family defaults to 0
#' @param Sib_status indicates the number of sib-cases defaults to 0
#' @keywords DistPlot
#' @export
#' @importFrom dplyr %>%
#' @examples
#' plosteriot(total_indiv = 1000, SNP = 1000, h = 0.5, Child=1, Mom=1, Dad=1, Nr_sib=0, Sib_status=0)

plosteriot <- function(total_indiv, SNP, h, Child=0, Mom=0, Dad=0, Nr_sib=0, Sib_status=0){
  res <- getMeanSD(total_indiv, SNP, h, Child, Mom, Dad, Nr_sib, Sib_status)
  if (is.list(res)){
    data.frame(x = c(-3, 3)) %>% ggplot2::ggplot(ggplot2::aes(x)) +
      ggplot2::geom_vline(xintercept = res$mean, linetype = "longdash", colour = "black", alpha = 1, size =1) +
      ggplot2::stat_function(fun = dnorm,
                             n = 5000,
                             position = "identity",
                             geom = "area",
                             fill = "blue",
                             alpha = 0.2,
                             args = list(mean = 0,
                                         sd = sqrt(h)),
                             size = 0.9, colour = "black") +
      ggplot2::stat_function(fun = dnorm,
                             position = "identity",
                             geom = "area",
                             fill = "red",
                             alpha = 0.2,
                             n = 5000,
                             args = list(mean = res$mean,
                                         sd = res$sd),
                             size = 0.9, colour = "black") +
      ggplot2::theme_bw() +
      ggplot2::labs(x = 'Genetic liability',
                    y = 'Density',
                    title = paste("Input of ",format(total_indiv,scientific = F)," individuals and ",format(SNP,scientific = F)," SNPs with ",h*100," heritability", sep="")) +
      ggplot2::scale_x_continuous(breaks = seq(-3, 3))
  }
  else if (res == -1){
    print("Please choose a different combination of cases, since the one you have entered is not present in data!")
  }
  else {
    print("No data exists! - Try another predefined or run our simulation")
  }
}
