% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moveRelevantFiles.R
\name{moveRelevantFiles}
\alias{moveRelevantFiles}
\title{moveRelevantFiles function}
\usage{
moveRelevantFiles(total_indiv, SNP, h, k)
}
\arguments{
\item{total_indiv}{Is the number of individuals for the data-file}

\item{SNP}{Is the number of SNPs for the data-file}

\item{h}{Is the heritability usually 0.5}

\item{k}{siginificance level - set to 0.05}
}
\description{
This function moves relevant files
}
\examples{
moveRelevantFiles(total_indiv = 1000, SNP = 1000, h = 0.5, k = 0.05)
}
\keyword{Bed}
\keyword{Make}
