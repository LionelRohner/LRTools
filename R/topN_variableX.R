#' Find top most variable row entries in a microarray data matrix
#'
#' This function can be used to extract the rows that have the greatest variance for further PCA analysis.
#'
#' @param mat A matrix with numerical values.
#' @param topN Upper bound of the number of most variable row entries. By default the 1000 most variable entries are considered.
#' @param plot Show or hide plot. Is set to "TRUE" by default.
#'
#' @return The input data.frame or matrix containing only the rows that have the highest variance in the data set.
#' @export
#'
#' @examples A <- data.frame(rnorm(1000,0.2,0.01),runif(1000,0.2,0.75),rnorm(1000,0.75,0.05),
#'                           rnorm(1000,0.2,0.01),runif(1000,0.2,0.6),rnorm(1000,0.6,0.05))
#'
#' rownames(A) <- paste("cg",sample(10000000:99999999,1000,replace = FALSE), sep ="")
#' colnames(A) <- rep(paste("Sample_",1:6,sep=""))
#'
#' A <- as.matrix(A)
#'
#' topN_variableX(mat = A,
#'                topN = 1000,
#'                plot = TRUE)
topN_variableX <- function(mat, topN = 1000, plot = T){

  # calculate row var and sort
  Var <- matrixStats::rowVars(mat)
  names(Var) <- rownames(mat)
  Var <- sort(Var, decreasing = T)

  # find top X
  VarN <- Var[c(1:topN)]

  if (plot == T){
    plot(VarN, main = paste("Top ", topN, "most variable rows in the matrix"))
  }

  mat_out <- mat[names(VarN),]
}
