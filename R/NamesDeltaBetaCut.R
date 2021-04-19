#' Extraction of CpG-probes IDs or indices based on delta beta-values
#'
#' The purpose of this function is to find the row names or indices of the CpG probe that have a certain absolute delta beta-value difference between two groups of choice. The function is restricted to comparison of two groups.
#'
#' @param pd A \emph{minfi}- or \emph{ChAMP}-type sample sheet.
#' @param feature The name of the group variable as it appears in the sample sheet.
#' @param beta A beta-value matrix, where rows are CpG probes and columns are samples.
#' @param group1_Name The string representation of the factor corresponding to group 1 in the group variable "feature".
#' @param group2_Name Same as above, but for the other group.
#' @param deltaBetaThreshold The desired delta beta-value threshold. The delta beta-value is the mean difference in beta-values between two groups. The value is converted to an absolute value. The default is 0.2.
#' @param returnIndeces By setting "returnIndices" to "TRUE" the function will return the indices of the CpG probes in the original beta-value matrix instead of the row names. Is set to "FALSE" by default.
#' @param sampleNames If the column containing the sample names in the pd is not named "Sample_Name" (ChAMP default naming convention), the column name has to be added separately. By default the function assumes the name to be "Sample_Name".
#'
#' @return A vector containing the row names (by default the CpG-probe ID) or indices of CpG-probes that have a user-defined delta beta difference between two groups of interest.
#' @export
#'
#' @examples # create some random matrix
#' A <- data.frame(runif(1000,0.2,1),runif(1000,0.3,1),runif(1000,0.5,0.99),
#'                 runif(1000,0.1,0.7),runif(1000,0.2,0.6),runif(1000,0,0.6))
#'
#' rownames(A) <- paste("cg",sample(10000000:99999999,1000,replace = FALSE), sep ="")
#' colnames(A) <- rep(paste("Sample_",1:6,sep=""))
#'
#' pd <- cbind.data.frame(colnames(A),c(rep("Group_1",3),rep("Group_2",3)))
#' colnames(pd) <- c("Sample_Names","Groups")
#'
#' inclusion <- NamesDeltaBetaCut(pd = pd,
#'                                feature = pd$Groups,
#'                                beta = A,
#'                                group1_Name = "Group_1",
#'                                group2_Name = "Group_2",
#'                                deltaBetaThreshold = 0.2,
#'                                returnIndeces = FALSE,
#'                                sampleNames = "Sample_Names")
#'
#' # Plot results
#' A_filtered <- A[rownames(A) %in% inclusion,]
#' Group1 <- paste("^",pd[which(pd$Groups == "Group_1"), "Sample_Names"],"$",
#'                 collapse = "|", sep = "")
#' Group2 <- paste("^",pd[which(pd$Groups == "Group_2"), "Sample_Names"],"$",
#'                 collapse = "|", sep = "")
#'
#' meanG1 <- rowMeans(A[,grep(Group1, colnames(A))])
#' meanG2 <- rowMeans(A[,grep(Group2, colnames(A))])
#'
#' hist(meanG1-meanG2, xlab = "delta beta-values",
#'      main = "All (grey) vs. filtered delta beta-values (red)")
#'
#' meanG1_filtered <- rowMeans(A_filtered[,grep(Group1, colnames(A_filtered))])
#' meanG2_filtered <- rowMeans(A_filtered[,grep(Group2, colnames(A_filtered))])
#'
#' hist(meanG1_filtered-meanG2_filtered, col = rgb(0.8,0.1,0.1,0.3),add = TRUE)

NamesDeltaBetaCut <- function(pd,
                              feature,
                              beta,
                              group1_Name,
                              group2_Name,
                              deltaBetaThreshold = 0.2,
                              returnIndeces = F,
                              sampleNames = "Sample_Name"){

  # Check whether sample names can be found
  `%notin%` <- Negate(`%in%`)
  if("Sample_Name" %notin% colnames(pd)){
    if(sampleNames %notin% colnames(pd)){
      message("Could not find variable \"Sample_Names\" in your pd. Please specify the exact variable name in argument \"sampleNames\". \nExiting...")
      return(NULL)
    }
  }

  # Mark delta-beta > |0.2| probes for filtration after conversion to M-values

  Group1 <- paste("^",pd[which(feature == group1_Name), sampleNames],"$",
                  collapse = "|", sep = "")

  Group2 <- paste("^",pd[which(feature == group2_Name), sampleNames],"$",
                  collapse = "|", sep = "")

  # Print Groups
  cat(group1_Name, ":", Group1, "\n")
  cat(group2_Name, ":", Group2, "\n")

  # Find mean difference
  Mean_Group1 <- rowMeans(beta[,grep(Group1,colnames(beta))])
  Mean_Group2 <- rowMeans(beta[,grep(Group2,colnames(beta))])


  # find name of the probes that are above delta-Beta threshold (absolute Value)
  if (returnIndeces == F){
    cut <- names(which((abs(Mean_Group1 - Mean_Group2) > deltaBetaThreshold) == T))
  } else{
    cut <- which((abs(Mean_Group1 - Mean_Group2) > deltaBetaThreshold) == T)

  }

  return(cut)
}
