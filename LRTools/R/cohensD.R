#' Cohen's D
#'
#' This function calculates Cohen's D, an effect size measure. If sample sizes are not equal use the unbiased modifier.
#'
#' @param meanA Mean of Group A
#' @param meanB Mean of Group B
#' @param sdA SD of Group A
#' @param sdB SD of Group B
#' @param nA Sample size of Group A
#' @param nB Sample size of Group B
#' @param unbiased If False the sample size will be assumed to be equal
#' @param verbose Extra printed information
#'
#' @return Cohen's D
#' @export
#'
#' @examples GroupA <- rnorm(1000,10, 2)
#' GroupB <- rnorm(1000,20, 4)
#'
#' meanGroupA <- mean(GroupA)
#' meanGroupB <- mean(GroupB)
#'
#' sigmaGroupA <- sd(GroupA)
#' sigmaGroupB <- sd(GroupB)
#'
#' cohensD(meanA = meanGroupA,
#'         meanB = meanGroupB,
#'         sdA = sigmaGroupA,
#'         sdB = sigmaGroupB,
#'         nA = 1000, nB = 1000,
#'         unbiased = TRUE,
#'         verbose = TRUE)
cohensD <- function(meanA,meanB,sdA,sdB,nA,nB, unbiased = F, verbose = F){


  if (verbose == T){
    spool <- sqrt((sdA^2+sdB^2)/2)
    print(paste("Pooled Variance : ", as.numeric(spool)),quote = F)
  }

  # unbiased has to be used if nA and nB are not equal
  if (unbiased == F){
    return((meanB-meanA)/(sqrt((sdA^2+sdB^2)/2)))
  } else {
    return((meanB-meanA)/(sqrt(((nA-1)*sdA^2 + (nB-1)*sdB^2)/(nA + nB - 2))))
  }
}


