#' Plotting PowerCalc Results
#'
#' Plotting function for powerCalc sample estimation. Generates a plot of the estimated sample size as a function of the CpG probes that have a power of 0.8.
#'
#' @param ... Any data.frames produced by powerCalc.
#' @param sampleSizeSteps Estimated sample sizes for which the number of CpG-probes with power of 0.8 are calculated. Represents the range that will be plotted. By default this is "c(5,10,15,20,25,30,35,40,45,50)".
#' @param returnCleanPlotObj Returns an empty ggplot2 object that can be modified. Is set to "FALSE" by default.
#' @param explicitNames The explicit names of the data.frames inputted in the dots argument. Maybe useful for plots shown in presentations. Assumes same order as the data.frames. Is not used by default.
#'
#' @return A plot or a ggplot object depending on returnCleanPlotObj choice.
#' @export
#'
#' @examples
#' # Standard Example
#'
#' A <- data.frame(runif(1000,0.5,1),runif(1000,0.4,0.8),runif(1000,0.6,0.99),
#'                 runif(1000,0.1,0.5),runif(1000,0.2,0.6),runif(1000,0,0.4))
#'
#' colnames(A) <- rep(paste("Sample_",1:6,sep=""))
#' pdA <- data.frame(c(rep("Group_1",3),rep("Group_2",3)))
#' colnames(pdA) <- "Groups"
#'
#' A <- as.matrix(A)
#'
#' powerCalcA <-powerCalc(betaMatrix = A,
#'                        type = "unpaired",
#'                        pdGroups = pdA$Groups,
#'                        nameGroups = c("Group_1","Group_2"),
#'                        cutOff = 0.1)
#'
#' B <- data.frame(runif(1000,0.8,1),runif(1000,0.8,0.9),runif(1000,0.7,0.9),
#'                 runif(1000,0.05,0.3),runif(1000,0.2,0.3),runif(1000,0,0.3))
#'
#' colnames(B) <- rep(paste("Sample_",1:6,sep=""))
#' pdB <- data.frame(c(rep("Group_1",3),rep("Group_2",3)))
#' colnames(pdB) <- "Groups"
#'
#' B <- as.matrix(B)
#'
#' powerCalcB <-powerCalc(betaMatrix = B,
#'                        type = "unpaired",
#'                        pdGroups = pdB$Groups,
#'                        nameGroups = c("Group_1","Group_2"),
#'                        cutOff = 0.1)
#'
#' powerCalcPlot(powerCalcA, powerCalcB)
#'
#' # Example with returnCleanPlotObj and explicitNames
#'
#' g <- powerCalcPlot(powerCalcA, powerCalcB,
#'                    returnCleanPlotObj = TRUE,
#'                    explicitNames = c("I Am","The Greatest"))
#' library(ggplot2)
#' g + geom_line(aes(color = Group)) + labs(title = "my Title") # + etc as you like
powerCalcPlot <- function(...,
                          sampleSizeSteps = c(5,10,15,20,25,30,35,40,45,50),
                          returnCleanPlotObj = FALSE,
                          explicitNames = NULL){

  # put all input dataframes in a list
  allPowerCalcObjects <- list(...)

  # names assigned to the data.frames
  if (is.null(explicitNames)){
    # get Arg names for group names (minus 1, as first arg is the function name)
    argNames <- as.list(sys.call())[-1]
  } else {
    argNames <- explicitNames
    if(length(argNames) != length(allPowerCalcObjects)){
      message("The number of names and data.frames do not match. \nExiting...")
      return(NULL)
    }
  }


  # INIT loop
  listResults <- list() # stores output data.frames
  cnt = 0 # used for naming and list appending

  for (powerCalcObj in allPowerCalcObjects){

    # counter used for correct naming
    cnt = cnt + 1

    # find number of probes that have power of 0.8 at different estimated N
    powerCalc_NvsProbes <-c()
    for (i in 1:length(sampleSizeSteps)){
      powerCalc_NvsProbes <- c(powerCalc_NvsProbes,
                               nrow(powerCalcObj[powerCalcObj$N <= sampleSizeSteps[i],]))
    }

    # add sample size steps
    powerCalc_NvsProbes <- cbind.data.frame(powerCalc_NvsProbes, sampleSizeSteps)

    # add analysis name as factor for plotting different dfs on same plot
    powerCalc_NvsProbes$G <- rep(factor(argNames[cnt]), nrow(powerCalc_NvsProbes))

    # add colnames for rbind in next step
    colnames(powerCalc_NvsProbes) <- c("Probes", "N", "Group")

    listResults[[cnt]] <- powerCalc_NvsProbes
  }


  # generate dataframe for ggplot
  df_4_ggplot <- data.frame()
  for (obj in listResults){
    df_4_ggplot <- rbind.data.frame(df_4_ggplot, obj)
  }

  # po-lo-ploto-ruuu
  if (returnCleanPlotObj == F){
    g <- ggplot2::ggplot(df_4_ggplot, ggplot2::aes(x = N, y = Probes, fill = Group)) +
      ggplot2::geom_line(ggplot2::aes(color = Group)) +
      ggplot2::geom_point(ggplot2::aes(color = Group), size = 2) +
      ggplot2::labs(title = "Power Analysis",
                    subtitle = "Sample size as a function of CpG probes with power 0.8") +
      ggplot2::xlab("Pairs of Sample") + ggplot2::ylab("CpG-Probes")

    print(g)
  } else {
    g <- ggplot2::ggplot(df_4_ggplot, ggplot2::aes(x = N, y = Probes, fill = Group))

    return(g)
  }
}
