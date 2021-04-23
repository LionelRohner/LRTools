#' Converts beta-values to M-values.
#' @export
#' @param beta numeric variable, vector or matrix
Beta_To_M <- function(beta){


  if (any(beta == 0) | any(beta == 1)){
    message("WARNING! \nYour matrix contains 0s or 1s, which turn to Inf during log-transform.")
    message("Please replace min and max values, such that they are not 0 or 1.")
    message("Exiting...")
    return(NULL)
  }


  return(log2(beta/(1-beta)))
}
