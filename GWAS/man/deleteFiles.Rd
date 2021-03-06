% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deleteFiles.R
\name{deleteFiles}
\alias{deleteFiles}
\title{deleteFiles function}
\usage{
deleteFiles(total_indiv, SNP, h, k)
}
\arguments{
\item{total_indiv}{Is the number of individuals for the data-file}

\item{SNP}{Is the number of SNPs for the data-file}

\item{h}{Is the heritability usually 0.5}

\item{k}{siginificance level - set to 0.05}
}
\description{
This function deletes relevant files
}
\examples{
deleteFiles(total_indiv = 100, SNP = 100, h = 0.5, k = 0.05)
}
\keyword{Bed}
\keyword{Make}
