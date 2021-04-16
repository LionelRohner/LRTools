#' powerCalc
#'
#'The purpose of this function is to evaluate the power and estimated sample size across all CpG probes present in a beta-value or M-value matrix. The calculation of power and sample size estimation are based on the t-test functions of the \emph{pwr} package and can only be performed between features that have two levels, e.g., primary PanNET versus metastasis. PowerCalc works for normal two-group designs or paired designs.
#'
#' @param betaMatrix A matrix of beta-values where the rows are the CpG probes and the columns are the samples.
#' @param pdGroups A \emph{minfi}- or \emph{ChAMP}-type sample sheet with group information, e.g. pd\$Sample\_Group.
#' @param nameGroups Name of the groups, e.g. c("Met", "Primary").
#' @param alpha Desired alpha-significance level.
#' @param power Desired power.
#' @param cutOff The power and sample size estimation is omitted if the effect size (Cohen's D) of a probe is less than the desired cutoff. Use this with caution as information about many CpG probes will be lost.
#' @param type Can be "paired" or "limma". Any other string will lead to the standard two.sample t-test procedure of \emph{pwr}.
#' @param limmaFit A \emph{LIMMA} LMfit object with one contrast can be added to calculate the effect size based on the posterior values for the residual variances obtained from the empirical Bayes procedure.
#' @param pairs Pair information of sample sheet, e.g. pd\$Pairs.
#' @param M Convert beta-values to M-values.
#'
#' @return A data.frame containing CpG probe-wise power, estimated sample size, and Cohen's D. Depending on the type other information are included (e.g. mean and SD of the difference between 2 groups when type = "paired).
#' @export
#'
#' @examples
#' A <- data.frame(runif(1000,0.5,1),runif(1000,0.4,0.8),runif(1000,0.6,0.99),
#'                 runif(1000,0.1,0.5),runif(1000,0.2,0.6),runif(1000,0,0.4))
#'
#' colnames(A) <- rep(paste("Sample_",1:6,sep=""))
#' pd <- data.frame(c(rep("Group_1",3),rep("Group_2",3)))
#' colnames(pd) <- "Groups"
#'
#' powerCalc_v3(betaMatrix = as.matrix(A),
#'              type = "unpaired",
#'              pdGroups = pd$Groups,
#'              nameGroups = c("Group_1","Group_2"),
#'              M = TRUE,
#'              cutOff = 0.1)
powerCalc_v3 <- function(betaMatrix,
                         pdGroups,
                         nameGroups,
                         alpha = 1e-6,
                         power = 0.8,
                         cutOff = 0,
                         type,
                         limmaFit = NULL,
                         pairs = NULL,
                         M = T){

  ### Step 0
  # Argument checking - Check if input are valid?

  # define notin
  `%notin%` <- Negate(`%in%`)

  message("")

  if (cutOff < 0.1){
    message("\nWARNING!\n")
    message("If the cutOff for effect size is too small the power analysis may crash, because")
    message("root finding is limited to the interval [1e-10, 1e+09]. To make sure what cutOff")
    message("is best suited for you, run powerCalc on a subset of the beta matrix and inspect")
    message("how the effect size (Cohens D) relates to sample size.\n")
    check <- readline(prompt="Press [y] to continue or any other key to exit: ")
    message("")

    if (check != "y"){
      message("\nExiting...\n")
      return(NULL)
    }
  }

  # Argument checking
  if (nlevels(pdGroups) > 2){
    message("Factors with more than 2 levels are not implemented. \nExiting...")
    return(NULL)
  }

  if (length(nameGroups) > 2 | length(nameGroups) == 1){
    message("Factors with more than 2 levels are not implemented. \nExiting...")
    return(NULL)
  }

  if (sum(nameGroups %notin% pdGroups) > 0){
    message("Group names do not fit sample sheet entries, review argument \"nameGroups\". \nExiting...")
    return(NULL)
  }


  if (type == "paired"){
    if(is.null(pairs)){
      message("Please add pair info of the sample sheet in argument \"pairs\".\nExiting...")
      return(NULL)
    }

    if (type == "limma"){
      if(is.null(limmaFit)){
        message("Please add a limma LMfit object with one contrast! \"pairs\".\nExiting...")
        return(NULL)
      }
    }

    if(length(unique(pairs)) != length(pairs)/2){
      message("Your samples pairs seem to be unbalanced, please check your sample sheet.\nExiting...")
      return(NULL)
    }
  }

  message("#######################################################################################")
  message("")
  message("||||||+   ||||||+  ||+    ||+ |||||||+ ||||||+      ||||||+  |||||+  ||+       ||||||+")
  message("||+--||+ ||+---||+ |||    ||| ||+----+ ||+--||+    ||+----+ ||+--||+ |||      ||+----+")
  message("||||||++ |||   ||| ||| |+ ||| |||||+   ||||||++    |||      |||||||| |||      |||     ")
  message("||+---+  |||   ||| ||||||+||| ||+--+   ||+--||+    |||      ||+--||| |||      |||     ")
  message("|||      +||||||++ +|||+|||++ |||||||+ |||  |||    +||||||+ |||  ||| |||||||+ +||||||+")
  message("+-+       +-----+   +--++--+  +------+ +-+  +-+     +-----+ +-+  +-+ +------+  +-----+")
  message("")
  message("#######################################################################################")
  message("")


  ### STEP 1
  # calculate means, SD of probes


  message("[Step 1]\n")
  if (type == "paired"){

    message("Paired T-Test was chosen. \n")
    message("Computing probe-wise mean and SD of the difference across pairs. \n")


    if (M == T){
      message("Converting beta- to M-values.\n")
      beta <- Beta_To_M(betaMatrix)
    }

    # add pairnames to colnames to guarantee pair-wise subtraction
    cNames <- paste(pairs,colnames(beta), sep = "_")

    colnames(beta) <- cNames

    # get groups
    M1 <- beta[,which(pdGroups == nameGroups[1])]
    M2 <- beta[,which(pdGroups == nameGroups[2])]

    # order samples accoding to pairs
    M1 <- M1[,order(colnames(M1))]
    M2 <- M2[,order(colnames(M2))]


    # REVIEW PAIR SUBTRACTION
    test <- cbind.data.frame(colnames(M2),colnames(M1))
    colnames(test) <- c(nameGroups[2], nameGroups[1])
    print(test)

    message("\nAre the sample pairs correct?\n")
    check <- readline(prompt="Press [y] to proceed with subtraction or any other key to exit: ")
    message("")
    if (check != "y"){
      message("Exiting...\n")
      return(NULL)
    }

    # Calculate difference of pairs (for paired analysis)
    D <- M2 - M1

    # mean of difference
    meanD <- rowMeans(D)

    # sd of difference
    sigmaD <- sqrt(matrixStats::rowVars(D))

    message("")
    message("[Step 2]\n")
    message("Computing effect size (Cohens D).\n")

    # calculate CohensD
    CohensD <- cohensD_Paired(meanD = meanD,
                              sigmaD = sigmaD)

    df_4_power <- cbind.data.frame(meanD,sigmaD,CohensD)


  } else if (type == "limma"){

    message("Limma estimates chosen. \n")

    if (M == T){
      message("Converting beta- to M-values.\n")
      beta <- Beta_To_M(betaMatrix)
    }

    CohensD <- (limmaFit$coefficients / sqrt(limmaFit$s2.post)) / sqrt(2) # paired ANALYSIS

    CohensD <- CohensD[rownames(CohensD) %in% rownames(beta),]

    df_4_power <- cbind.data.frame(CohensD)


  } else{
    message("Two-Sample T-Test was chosen. \n")
    message("Computing Mean and SD of probes for each group. \n")

    if (M == T){
      message("Converting beta- to M-values.\n")
      beta <- Beta_To_M(betaMatrix)
    }

    # mean
    M1 <- rowMeans(beta[,which(pdGroups == nameGroups[1])])
    M2 <- rowMeans(beta[,which(pdGroups == nameGroups[2])])

    SD1 <- sqrt(matrixStats::rowVars(beta[,which(pdGroups == nameGroups[1])]))
    SD2 <- sqrt(matrixStats::rowVars(beta[,which(pdGroups == nameGroups[2])]))

    # N
    n1 = length(which(pdGroups == nameGroups[1]))
    n2 = length(which(pdGroups == nameGroups[2]))


    # create DF
    message("")
    message("[Step 2]\n")
    message("Computing effect size (Cohens D).\n")
    message("Mean difference : Mean(Group2) - Mean(Group1)")

    df_4_power <- cbind.data.frame(M1,M2,SD1,SD2)
    df_4_power$meanDiff <- (df_4_power$M2 - df_4_power$M1)

    # compute cohens D
    df_4_power$CohensD <- cohensD(meanA = df_4_power$M1,
                                  meanB = df_4_power$M2,
                                  sdA = df_4_power$SD1,
                                  sdB = df_4_power$SD2,
                                  nA = n1, nB = n2,
                                  unbiased = T)
  }


  ### Step 2 - Calculate probe-wise estimated sample size and power


  # loop takes very long
  powerAnalysis <- data.frame()


  # remove up to cutOff
  df_4_power <- df_4_power[which(abs(df_4_power$CohensD) > cutOff),]

  # UGLY WORK AROUND
  if (type == "limma"){
    df_4_power <- as.data.frame(df_4_power)
    colnames(df_4_power) <- c("CohensD")
  }


  loop_lim <- nrow(df_4_power)

  message("")
  message("[Step 3]\n")
  message("Cohens D Cutoff is set to ",cutOff,". Hence, ",round((1-nrow(df_4_power)/nrow(beta))*100, digits = 3),"% of the probes will be removed form the analysis.\n")
  message("PowerCalc in progress... \n")
  message("Loop might get slower towards the end, algo is not perfectly vectorized. \n\n")


  # change mods depending on analysis
  if (type == "paired"){
    typeMod <- "paired"
    n <- length(unique(pairs))

  } else if (type == "limma"){
    typeMod <- "paired"
    n <- length(unique(pairs))

  } else {
    typeMod <- "two.sample"
    n <- length(pdGroups)
  }

  for (i in 1:loop_lim){
    CohensD = df_4_power$CohensD[i]

    res <- c(pwr::pwr.t.test(d = CohensD,
                             type = typeMod, alternative = ("two.sided"),
                             power = 0.8,
                             sig.level = alpha)$n,

             pwr::pwr.t.test(d = CohensD,
                             type = typeMod,
                             n = n,
                             sig.level = alpha)$power,
                             CohensD)


    # populate output
    powerAnalysis = rbind(powerAnalysis, res)

    # progress
    cat(sprintf("\r%.3f%%", (i / (loop_lim) * 100)))
  }
  message("\n")
  message("[Step 4]\n")
  message("Formatting output.\n")

  colnames(powerAnalysis) = c("N","Power","Cohens D")

  pA <- cbind(df_4_power, powerAnalysis)

  message("\nReturning Dataframe...")
  return(pA)
}
