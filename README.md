# LRTools
## Introduction ##
This is a small collection of the most useful functions I wrote during my master project. The package consists of functions that are necessary for the bioinformatics pipeline that I used to analyze Illumina HM450 and EPIC data. The package also includes other useful functions that can be helpful when analyzing DNA methylation data from these platforms.

## Installation & Dependecies ##
To install this package use : ```devtools::install_github("LionelRohner/LRTools")``` <br />
Following packages are required to use all functions of this pacakge: ```matrixStats, pwr, ChAMP, ChAMPdata, akmedoids, ggplot2```

## Functions ##

### Beta_To_M ###

Conversion of beta-values to M-values.


### champ.load_extended ###

Modified champ.load function from the \emph{ChAMP} package.


### checkBMIQ ###

Plot the beta distribution of Infinium I and II probes before and after BMIQ-correction.


### cohensD_Paired ###

Calculate Cohen's D effect size for paired data. 


### cohensD ###

Calculate standard Cohen's D effect size.


### exVarPlot ###

Visualize explained variance of PC-1 and PC-2 of a matrix contaning only the most variable rows.


### NamesDeltaBetaCut ###

Extraction of CpG-probes IDs or indices based on delta beta-values

### powerCalc ###

Power and sample size estimation.



### powerCalcPlot ###

Plotting PowerCalc Results.


### topN_variableX ###

Find top most variable row entries in a microarray data matrix.
