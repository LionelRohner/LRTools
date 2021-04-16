#' Visualize explained variance of PC-1 and PC-2 of a matrix contaning only the most variable rows.
#'
#' To be honest, I am not entirely sure how informative this function is. I used it to see how much the explained variance of the PCs varies in the lower tail in order to choose a suitable number of the most variable CpG-probes. This number can be quite arbitrary, and if this number is in a region of the distribution that fluctuates a lot, it may be wiser to reconsider the number of most variable rows. Explained variance is calculated as the ratio of eigenvalues of PC-1 or PC-2 divided by the sum of all eigenvalues as proposed here: https://stats.stackexchange.com/questions/254592/calculating-pca-variance-explained/254598.
#'
#' @param mat A matrix of numeric values, where rows are probes, genes, metabolites, etc. and the columns are samples.
#' @param lowerLim The lower limit of the most variable row entries.
#' @param upperLim The upper limit of the most variable row entries.
#' @param stepSize The step size. By default the step size is 10.
#' @param onlyPC1 Only consider PC-1. Is "FALSE" by default.
#'
#' @return A plot and a data.frame with the explained variance of PC-1 or PC-1 + PC-2 at each step.
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
#' exVarPlot(mat=as.matrix(A),
#'           lowerLim = 10,
#'           upperLim = 1000)
exVarPlot <- function(mat, lowerLim, upperLim, stepSize = 10, onlyPC1 = F){

  # TO DO: remove dependency
  if (requireNamespace("akmedoids") == F){
    message("Please install \"akmedoids\" first. Exiting...\n")
    return(NULL)
  }

  # empy objects for results
  index <- c()
  varEx <- c()

  # get top N most variable probes
  topN <- t(topN_variableX((mat), topN = upperLim, plot = F))

  message("Calculating explained variance of PC1 and PC2 for different N most variable probes/genes.")
  for (i in seq(lowerLim, upperLim, stepSize)){

    pca_res <- prcomp(topN[,1:i], scale. = TRUE)

    # calculate explained variance in PC1 and PC2 and add them up.
    # Doing so produces less fluctuation in the data compared to using PC1 only
    eigs <- pca_res$sdev^2

    if (onlyPC1 == F){
      varEx <- append(varEx,(eigs[1] / sum(eigs) + eigs[2] / sum(eigs)))
      index <- append(index, i)
    } else {
      varEx <- append(varEx,(eigs[1] / sum(eigs)))
      index <- append(index, i)
    }

    # progress
    cat(sprintf("\r%.3f%%", (i / (upperLim) * 100)))
  }
  message("\n")
  # save results
  Res <- cbind.data.frame(index,varEx)

  message("Performing LOESS curve fitting.\n")

  # fit loess to smooth the curve
  loessCurve <- loess(varEx ~ index,data=Res)

  j <- order(Res$index)
  elB <- cbind.data.frame(Res$index[j],loessCurve$fitted[j])

  # compute elbow point
  elbowPnt <- akmedoids::elbow_point(elB[,1],elB[,2])

  # plot
  if (onlyPC1 == F){
    plot(varEx ~ index, data=Res,pch=19,cex=0.1,
         main = "Explained Variance of the Sum of PC1 and PC2 of the Most Variable Probes",
         xlab = "Index", ylab = "Var. Explained [PC1 + PC2]")
  } else {
    plot(varEx ~ index, data=Res,pch=19,cex=0.1,
         main = "Explained Variance of the PC1 of the Most Variable Probes",
         xlab = "Index", ylab = "Var. Explained [PC1]")
  }

  abline(v = elbowPnt$x, lty = 2, col = rgb(0.5,0.5,0.5,0.5))
  abline(h = elbowPnt$y, lty = 2, col = rgb(0.5,0.5,0.5,0.5))
  points(elbowPnt$x,elbowPnt$y, col = rgb(0.5,0.5,0.5,0.5))
  legend("topright", inset = 0.02,
         legend=c("LOESS",paste("Elbow Point: ",round(elbowPnt$x))),
         lty = 1:2,col=c(rgb(0.9,0.3,0.1,0.7), rgb(0.5,0.5,0.5,0.5)))

  lines(Res$index[j],loessCurve$fitted[j],col=rgb(0.9,0.1,0.1,0.7),lwd=2)

  message("Returning results...")
  return(Res)
}
