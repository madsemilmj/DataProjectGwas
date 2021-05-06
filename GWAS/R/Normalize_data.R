#' Normalize_data function
#'
#' This function assigns phenotypes based on GWAX.
#' @param df Is a matrix/data-frame that we want to Normalize (number_of_indivs X SNPs)
#' @param MAF Defines the Minor Allele Frequancy, and is a vector of probabilities from unif(0.01,0.49) one for each SNP
#' @keywords Normalize
#' @export
#' @examples
#' Normalize_data()


Normalize_data <- function(df,MAF){
  df_normal <- sweep(sweep(df, 2, 2*MAF, "-"), 2, sqrt(2*MAF*(1-MAF)), "/")
  return(df_normal)
}
