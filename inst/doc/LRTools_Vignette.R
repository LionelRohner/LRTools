## ---- include = FALSE-----------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = F---------------------------------------------------------------------------------------------------------------------
library(LRTools)

## -------------------------------------------------------------------------------------------------------------------------------------------
vecBeta <- runif(100,0,1)
vecM <- Beta_To_M(vecBeta)
head(vecM)

## ---- eval=FALSE----------------------------------------------------------------------------------------------------------------------------
#  myLoad_450K <- champ.load_extended(method = "minfi",
#                                   arraytype = "450K",
#                                   sampleSheet.csv ="^someSampleSheet.csv$",
#                                   preproc = "Noob",
#                                   dyeMethod = "single")

## ---- eval=FALSE----------------------------------------------------------------------------------------------------------------------------
#  checkBMIQ(betaBefore = beta,
#                      betaAfter = beta_BMIQcorrected,
#                      array = "450K")

## -------------------------------------------------------------------------------------------------------------------------------------------
before <- rnorm(1000, 10, 2)
after <- rnorm(1000,20, 4)
D <- before - after
meanD <- mean(D)
sigmaD <- sd(D)
cohensD_Paired(meanD = meanD,
               sigmaD = sigmaD)


## -------------------------------------------------------------------------------------------------------------------------------------------
GroupA <- rnorm(1000,10, 2)
GroupB <- rnorm(1000,20, 4)

meanGroupA <- mean(GroupA)
meanGroupB <- mean(GroupB)

sigmaGroupA <- sd(GroupA)
sigmaGroupB <- sd(GroupB)

cohensD(meanA = meanGroupA,
        meanB = meanGroupB,
        sdA = sigmaGroupA,
        sdB = sigmaGroupB,
        nA = 1000, nB = 1000,
        unbiased = TRUE,
        verbose = TRUE)


## ---- results='hide', message=FALSE,fig.width = 8, fig.height = 5---------------------------------------------------------------------------
A <- data.frame(rnorm(1000,0.2,0.01),runif(1000,0.2,0.75),rnorm(1000,0.75,0.05),
                          rnorm(1000,0.2,0.01),runif(1000,0.2,0.6),rnorm(1000,0.6,0.05))

rownames(A) <- paste("cg",sample(10000000:99999999,1000,replace = FALSE), sep ="")
colnames(A) <- rep(paste("Sample_",1:6,sep=""))

A <- as.matrix(A)

X <- exVarPlot(mat=A,
               lowerLim = 10,
               upperLim = 1000)


## -------------------------------------------------------------------------------------------------------------------------------------------
head(X)

## ---- fig.width = 8, fig.height = 5---------------------------------------------------------------------------------------------------------
# create some random matrix
A <- data.frame(runif(1000,0.2,1),runif(1000,0.3,1),runif(1000,0.5,0.99),
                runif(1000,0.1,0.7),runif(1000,0.2,0.6),runif(1000,0,0.6))

rownames(A) <- paste("cg",sample(10000000:99999999,1000,replace = FALSE), sep ="")
colnames(A) <- rep(paste("Sample_",1:6,sep=""))

pd <- cbind.data.frame(colnames(A),c(rep("Group_1",3),rep("Group_2",3)))
colnames(pd) <- c("Sample_Names","Groups")

inclusion <- NamesDeltaBetaCut(pd = pd,
                               feature = pd$Groups,
                               beta = A,
                               group1_Name = "Group_1",
                               group2_Name = "Group_2",
                               deltaBetaThreshold = 0.2,
                               returnIndeces = FALSE,
                               sampleNames = "Sample_Names")

# Plot results
A_filtered <- A[rownames(A) %in% inclusion,]
Group1 <- paste("^",pd[which(pd$Groups == "Group_1"), "Sample_Names"],"$",
                collapse = "|", sep = "")
Group2 <- paste("^",pd[which(pd$Groups == "Group_2"), "Sample_Names"],"$",
                collapse = "|", sep = "")

meanG1 <- rowMeans(A[,grep(Group1, colnames(A))])
meanG2 <- rowMeans(A[,grep(Group2, colnames(A))])

hist(meanG1-meanG2, xlab = "delta beta-values",
     main = "All (grey) vs. filtered delta beta-values (red)")

meanG1_filtered <- rowMeans(A_filtered[,grep(Group1, colnames(A_filtered))])
meanG2_filtered <- rowMeans(A_filtered[,grep(Group2, colnames(A_filtered))])

hist(meanG1_filtered-meanG2_filtered, col = rgb(0.8,0.1,0.1,0.3),add = TRUE)

## ---- results=F, message=F, fig.width = 8, fig.height = 5-----------------------------------------------------------------------------------
A <- data.frame(runif(1000,0.5,1),runif(1000,0.4,0.8),runif(1000,0.6,0.99),
                runif(1000,0.1,0.5),runif(1000,0.2,0.6),runif(1000,0,0.4))

colnames(A) <- rep(paste("Sample_",1:6,sep=""))
pd <- data.frame(c(rep("Group_1",3),rep("Group_2",3)))
colnames(pd) <- "Groups"

A <- as.matrix(A)

PC <- powerCalc(betaMatrix = A,
                type = "unpaired",
                pdGroups = pd$Groups,
                nameGroups = c("Group_1","Group_2"),
                M = TRUE,
                cutOff = 0.1)

## -------------------------------------------------------------------------------------------------------------------------------------------
head(PC)

## ---- results=F, message=F, fig.width = 8, fig.height = 5-----------------------------------------------------------------------------------
# Standard Example

A <- data.frame(runif(1000,0.5,1),runif(1000,0.4,0.8),runif(1000,0.6,0.99),
                runif(1000,0.1,0.5),runif(1000,0.2,0.6),runif(1000,0,0.4))

colnames(A) <- rep(paste("Sample_",1:6,sep=""))
pdA <- data.frame(c(rep("Group_1",3),rep("Group_2",3)))
colnames(pdA) <- "Groups"

A <- as.matrix(A)

powerCalcA <-powerCalc(betaMatrix = A,
                       type = "unpaired",
                       pdGroups = pdA$Groups,
                       nameGroups = c("Group_1","Group_2"),
                       cutOff = 0.1)

B <- data.frame(runif(1000,0.8,1),runif(1000,0.8,0.9),runif(1000,0.7,0.9),
                runif(1000,0.05,0.3),runif(1000,0.2,0.3),runif(1000,0,0.3))

colnames(B) <- rep(paste("Sample_",1:6,sep=""))
pdB <- data.frame(c(rep("Group_1",3),rep("Group_2",3)))
colnames(pdB) <- "Groups"

B <- as.matrix(B)

powerCalcB <-powerCalc(betaMatrix = B,
                       type = "unpaired",
                       pdGroups = pdB$Groups,
                       nameGroups = c("Group_1","Group_2"),
                       cutOff = 0.1)

powerCalcPlot(powerCalcA, powerCalcB)

# Example with returnCleanPlotObj and explicitNames

g <- powerCalcPlot(powerCalcA, powerCalcB,
                   returnCleanPlotObj = TRUE,
                   explicitNames = c("Some","Thing"))
library(ggplot2)
g + geom_line(aes(color = Group)) + labs(title = "my Title") # + etc as you like

## ---- results=T, message=F, fig.width = 8, fig.height = 5-----------------------------------------------------------------------------------
A <- data.frame(rnorm(1000,0.2,0.01),runif(1000,0.2,0.75),rnorm(1000,0.75,0.05),
                          rnorm(1000,0.2,0.01),runif(1000,0.2,0.6),rnorm(1000,0.6,0.05))

rownames(A) <- paste("cg",sample(10000000:99999999,1000,replace = FALSE), sep ="")
colnames(A) <- rep(paste("Sample_",1:6,sep=""))

A <- as.matrix(A)

top1K <- topN_variableX(mat = A,
                        topN = 1000,
                        plot = TRUE)

head(top1K)

