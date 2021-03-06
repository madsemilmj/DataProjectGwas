% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MasterFunc.R
\name{MasterFunc}
\alias{MasterFunc}
\title{Master Func function}
\usage{
MasterFunc(total_indiv, indiv_chunk, SNP, h, c, k)
}
\arguments{
\item{total_indiv}{The total number of individuals to simulate}

\item{indiv_chunk}{The number of individuals to simulate in each iteration (should not be above 5000)}

\item{SNP}{The number of SNPs to simulate}

\item{h}{the h^2 - meaning the .. ususally 0.5}

\item{c}{The number of causal SNPs in the simulaiton usually 1/10000*SNP}

\item{k}{Significance level - usually 0.05}
}
\description{
This function Runs the simulation and assigns differet phenotypes (GWAX, LTFH) and generete relevant files used in furhter analysis
}
\examples{
MasterFunc(1000,100,1000,0.5,100,0.05)
}
\keyword{Master}
