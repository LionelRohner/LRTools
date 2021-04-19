#' Cohen's D for paired data
#'
#' This function calculates Cohen's D for paired data. The formula was taken from the \emph{rstatix} package.
#'
#' @param meanD Mean of the differences between measurements
#' @param sigmaD SD of the differences between measurements
#'
#' @return Cohen's D for paired data. Adapted from \emph{rstatix}
#' @export
#'
#' @examples before <- rnorm(1000, 10, 2)
#' after <- rnorm(1000,20, 4)
#' D <- before - after
#' meanD <- mean(D)
#' sigmaD <- sd(D)
#' cohensD_Paired(meanD = meanD,
#'                sigmaD = sigmaD)
cohensD_Paired <- function(meanD, sigmaD){
  return(meanD/sigmaD)
}
