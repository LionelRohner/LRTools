#' Conversion of beta-values to M-values
#'
#' This function performs a logit transform (\eqn{M=log2(beta)-log2(beta-1)}) on beta-values as single values or in vector, matrix, or data.frame form. M-values are not defined for beta-values equal to 0 or 1. A warning message will be displayed in case of presence of 0 and 1 in the matrix.
#'
#' @param beta A numeric variable, vector or matrix corresponding to a beta-value. Has to be between 0 and 1.
#'
#' @return logit transform of the input beta-values.
#' @export
#'
#' @examples
#' vecBeta <- runif(100,0,1)
#' Beta_To_M(vecBeta)
Beta_To_M <- function(beta){

  # prevent NaNs if 0 or 1 values are present
    if (any(beta == 0) | any(beta == 1)){
    message("WARNING! \nYour matrix contains 0s or 1s, which turn to Inf during log-transform.")
    message("Please replace min and max values, such that they are not 0 or 1.")
    message("Exiting...")
    return(NULL)
  }
  return(log2(beta/(1-beta)))
}
