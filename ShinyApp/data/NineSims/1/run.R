library(devtools)
library(tidyverse)
library(data.table)
library(stringr)
library(future.apply)
library(flock)
install_github("madsemilmj/DataProjectGwas/GWAS")

library(GWAS)

total_indiv <- 100000
indiv_chunk <- 100
SNP <- 100000
h <- 0.5
c <- SNP/100
k <- 0.05


MasterFunc(total_indiv = total_indiv, indiv_chunk = indiv_chunk, SNP = SNP, h = h, c = c, k = k)
