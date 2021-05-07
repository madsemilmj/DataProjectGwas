library(data.table)
library(tidyverse)
#' Plot distributions
#'
#' This function plots posterior distribution and prior
#' @param total_indiv Number of indicviduls in the dataset
#' @param SNP Is the number of SNPs in the data
#' @param h Is the heritability usually 0.5
#' @param Child Is a binary for indicating case/notcase for child
#' @param Mom Is a binary for indicating case/notcase for Mom
#' @param Dad Is a binary for indicating case/notcase for Dad
#' @param Nr-Sib indiacates the number of siblings for the family defaults to 0
#' @param Sib_status indicates the number of sib-cases defaults to 0
#' @keywords DistPlot
#' @export
#' @examples
#' plosteriot(total_indiv, SNP, h, Child, Mom, Dad, Nr_sib=0, Sib_status=0)


plosteriot <- function(total_indiv, SNP, h, Child, Mom, Dad, Nr_sib=0, Sib_status=0){
  res <- getMeanSD(total_indiv, SNP, h, Child, Mom, Dad, Nr_sib, Sib_status)
  if (res != -1){
    data.frame(x = c(-3, 3)) %>% ggplot(aes(x)) +
      geom_vline(xintercept = res$mean, linetype = "longdash", colour = "black", alpha = 1, size =1) +
      stat_function(fun = dnorm,
                    n = 5000,
                    position = "identity",
                    geom = "area",
                    fill = "blue",
                    alpha = 0.2,
                    args = list(mean = 0,
                                sd = 1),
                    size = 0.9, colour = "black") +
      stat_function(fun = dnorm,
                    position = "identity",
                    geom = "area",
                    fill = "red",
                    alpha = 0.2,
                    n = 5000,
                    args = list(mean = res$mean,
                                sd = res$sd),
                    size = 0.9, colour = "black") +
      theme_bw() +
      labs(x = 'Genetic liability',
          y = 'Density',
          title = 'Case with...') +
      scale_x_continuous(breaks = seq(-3, 3))
  }
  else{
    print("Please choose a different combination of cases, since the one you have entered is not present in data!")
  }
}
