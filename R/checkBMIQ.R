#' Plot the beta distribution of Infinium I and II probes before and after BMIQ-correction.
#'
#' This function can be used to evaluate the beta-value distribution of Infinium I and II probes before and after probe bias correction with BMIQ. Although implemented to check BMIQ performance, this function can be used to evaluate Infinium I and II probe distribution at any step of the pipeline and is also not limited to BMIQ-corrected beta-value matrices. By default, the function represents the mean density of Infinium I and II probes across all samples. It may be advisable to check for updates to \emph{ChAMP} and \emph{ChAMPdata} as new EPIC manifests are released from time to time.
#'
#' @param betaBefore The beta-value matrix before BMIQ correction. The rows should be CpG probes and the columns samples.
#' @param betaAfter The beta-value matrix after BMIQ correction. The rows should be CpG probes and the columns samples.
#' @param array Selection of the array type. For EPIC data use "EPIC" and everything else will be interpreted as HM450 data. By default HM450 data is assumed.
#' @param sampleChoice To check individual samples, enter a sample name as it appears in the column names. Multiple samples can be entered with regular expressions, e.g. ""aPxxx|aPyy"".
#'
#' @return No return. A plot should appear in the \emph{Plots} pane.
#' @export
#'
#' @examples \dontrun{checkBMIQ(betaBefore = beta,
#'                     betaAfter = beta_BMIQcorrected,
#'                     array = "450K")}
checkBMIQ <- function(betaBefore,
                      betaAfter,
                      array = "450K",
                      sampleChoice = NULL){

  # graphical params
  par(mfrow = c(1,2))

  # check packages
  if (requireNamespace("ChAMPdata") == F){
    message("Please install \"ChAMPdata\" first. Exiting...\n")
    return(NULL)
  }

  # check args
  if(!is.null(sampleChoice)){
    if (length(grep(sampleChoice,colnames(betaBefore))) == 0){
      message("Sample not found. Exiting...\n")
      par(mfrow = c(1,1))
      return(NULL)}
  }

  # check array type
  if (array == "EPIC")
    data(probeInfoALL.epic.lv)
  else data(probeInfoALL.lv)
  ProbeDesign <- as.numeric(lapply(probeInfoALL.lv,
                                   function(x) x)$Design[match(rownames(betaBefore),
                                                               probeInfoALL.lv$probeID)])
  # get indeces of type I and II probes in uncorrected matrix
  betaP1 <- betaBefore[which(ProbeDesign == 1),]
  betaP2 <- betaBefore[which(ProbeDesign == 2),]

  # plot before
  if (is.null(sampleChoice)){
    plot(density(rowMeans(betaP1)),xlab = "Beta",main = "Before BMIQ",
         cex.lab = 1.3, lwd = 1.5)
    lines(density(rowMeans(betaP2)), col = "firebrick1", lwd = 1.5)
    legend("topright", inset = 0.02, lty = 1,
           legend=c("Infinium I", "Infinium II"),
           col = c("black","firebrick1"))
  } else {
    plot(density(betaP1[,grep(sampleChoice,colnames(betaP1))]),xlab = "Beta",
         main = paste("Before BMIQ ","(",sampleChoice,")", sep = ""),
         cex.lab = 1.3, lwd = 1.5)
    lines(density(betaP2[,grep(sampleChoice,colnames(betaP2))]),
          col = "firebrick1", lwd = 1.5)
    legend("topright", inset = 0.02, lty = 1,
           legend=c("Infinium I", "Infinium II"),
           col = c("black","firebrick1"))
  }

  # get indices of type I and II probes in BMIQ-corrected matrix
  betaP1_after <- betaAfter[which(ProbeDesign == 1),]
  betaP2_after <- betaAfter[which(ProbeDesign == 2),]


  # plot after
  if (is.null(sampleChoice)){
    plot(density(rowMeans(betaP1_after)),xlab = "Beta",main = "After BMIQ",
         cex.lab = 1.3, lwd = 1.5)
    lines(density(rowMeans(betaP2_after)), col = "firebrick1", lwd = 1.5)
    legend("topright", inset = 0.02, lty = 1,
           legend=c("Infinium I", "Infinium II"),
           col = c("black","firebrick1"))
  } else {
    plot(density(betaP1_after[,grep(sampleChoice,colnames(betaP1_after))]),
         xlab = "Beta",
         main = paste("After BMIQ ","(",sampleChoice,")", sep = ""),
         cex.lab = 1.3, lwd = 1.5)
    lines(density(betaP2_after[,grep(sampleChoice,colnames(betaP2_after))]),
          col = "firebrick1", lwd = 1.5)
    legend("topright", inset = 0.02, lty = 1,
           legend=c("Infinium I", "Infinium II"),
           col = c("black","firebrick1"))
  }

  # Reminder
  message("ChAMPdata version ",package.version("ChAMPdata"),
          " was used. It is important to note that new versions ChAMP data may exist.")


  # undo split screen
  par(mfrow = c(1,1))
}
