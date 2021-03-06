% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topN_variableX.R
\name{topN_variableX}
\alias{topN_variableX}
\title{Find top most variable row entries in a microarray data matrix}
\usage{
topN_variableX(mat, topN = 1000, plot = T)
}
\arguments{
\item{mat}{A matrix with numerical values.}

\item{topN}{Upper bound of the number of most variable row entries. By default the 1000 most variable entries are considered.}

\item{plot}{Show or hide plot. Is set to "TRUE" by default.}
}
\value{
The input data.frame or matrix containing only the rows that have the highest variance in the data set.
}
\description{
This function can be used to extract the rows that have the greatest variance for further PCA analysis.
}
\examples{
A <- data.frame(rnorm(1000,0.2,0.01),runif(1000,0.2,0.75),rnorm(1000,0.75,0.05),
                          rnorm(1000,0.2,0.01),runif(1000,0.2,0.6),rnorm(1000,0.6,0.05))

rownames(A) <- paste("cg",sample(10000000:99999999,1000,replace = FALSE), sep ="")
colnames(A) <- rep(paste("Sample_",1:6,sep=""))

A <- as.matrix(A)

topN_variableX(mat = A,
               topN = 1000,
               plot = TRUE)
}
