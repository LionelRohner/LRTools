#' Modified champ.load function from the \emph{ChAMP} package
#'
#' This is a modification of the function "champ.load" from the \emph{ChAMP} package, which allows to use the Noob background correction. The original function from \emph{ChAMP} loads data from IDAT files to calculate intensity and produce quality control images. Furthermore, a new argument has been added that allows to load a specific sample sheet using a regular expression. This avoids the errors of the normal loading function when multiple csv files are found in the working directory. Only the new arguments are presented here. The other arguments can be found in the \emph{ChAMP} manual here: https://rdrr.io/bioc/ChAMP/man/champ.load.html.
#' The modification are minor, see line 94 and 104 to 111, everything else is from the \emph{ChAMP} package created by Yuan Tian (cre,aut), Tiffany Morris (ctb), Lee Stirling (ctb), Andrew Feber (ctb), Andrew Teschendorff (ctb), Ankur Chakravarthy (ctb). \emph{ChAMP} is licensed under GPL-3.
#'
#' @param preproc Choice of the preprocessing method of the rgSet. To use Noob dye bias and background correction use "Noob". The function is available in the \emph{minfi} package. Using "Raw" is the standard \emph{ChAMP} method and does not include background correction.
#' @param dyeMethod To date single sample procedure of Noob ("single") is the only choice here.
#' @param sampleSheet.csv An optional argument, that allows to add the name of a specific sample sheet using regular expressions.
#' @param directory See ChAMP manual.
#' @param method See ChAMP manual.
#' @param methValue See ChAMP manual.
#' @param autoimpute See ChAMP manual.
#' @param filterDetP See ChAMP manual.
#' @param ProbeCutoff See ChAMP manual.
#' @param SampleCutoff See ChAMP manual.
#' @param detPcut See ChAMP manual.
#' @param filterBeads See ChAMP manual.
#' @param beadCutoff See ChAMP manual.
#' @param filterNoCG See ChAMP manual.
#' @param filterSNPs See ChAMP manual.
#' @param population See ChAMP manual.
#' @param filterMultiHit See ChAMP manual.
#' @param filterXY See ChAMP manual.
#' @param force See ChAMP manual.
#' @param arraytype See ChAMP manual.
#'
#' @return A list including the rgSet, mSet, beta-value matrix, detection p-value matrix, and a matrix of intensity values. See \emph{ChAMP} manual for more detail.
#' @export
#'
#' @examples \dontrun{myLoad_450K <- champ.load_extended(method = "minfi",
#'                                  arraytype = "450K",
#'                                  sampleSheet.csv ="^someSampleSheet.csv$",
#'                                  preproc = "Noob",
#'                                  dyeMethod = "single")}
champ.load_extended <-function (directory = getwd(),
                                method = "ChAMP",
                                methValue = "B",
                                autoimpute = TRUE,
                                filterDetP = TRUE,
                                ProbeCutoff = 0,
                                SampleCutoff = 0.1,
                                detPcut = 0.01,
                                filterBeads = TRUE,
                                beadCutoff = 0.05,
                                filterNoCG = TRUE,
                                filterSNPs = TRUE,
                                population = NULL,
                                filterMultiHit = TRUE,
                                filterXY = TRUE,
                                force = FALSE,
                                arraytype = "450K",
                                sampleSheet.csv = "csv$",
                                preproc = "Noob",
                                dyeMethod = "single")
{
  message("[===================================]")
  message("[<<<< ChAMP.LOAD_EXTENDED START >>>>]")
  message("-------------------------------------")
  mybeadcount <- function(x) {
    nb <- getNBeads(x)
    typeIadd <- getProbeInfo(x, type = "I")
    typeImatchA <- match(typeIadd$AddressA, rownames(nb))
    typeImatchB <- match(typeIadd$AddressB, rownames(nb))
    typeIIadd <- getProbeInfo(x, type = "II")
    typeIImatch <- match(typeIIadd$Address, rownames(nb))
    nbcg <- nb
    locusNames <- getManifestInfo(x, "locusNames")
    bc_temp <- matrix(NA_real_, ncol = ncol(x), nrow = length(locusNames),
                      dimnames = list(locusNames, sampleNames(x)))
    TypeII.Name <- getProbeInfo(x, type = "II")$Name
    bc_temp[TypeII.Name, ] <- nbcg[getProbeInfo(x, type = "II")$AddressA,
    ]
    TypeI <- getProbeInfo(x, type = "I")
    bcB <- bc_temp
    bcA <- bc_temp
    bcB[TypeI$Name, ] <- nbcg[TypeI$AddressB, ]
    bcA[TypeI$Name, ] <- nbcg[TypeI$AddressA, ]
    bcB3 <- which(bcB < 3)
    bcA3 <- which(bcA < 3)
    bcA2 <- bcA
    bcB2 <- bcB
    bcA2[bcA3] <- NA
    bcA2[bcB3] <- NA
    bc <- data.frame(bcA2)
    bc
  }
  if (method == "minfi") {
    message("\n[ Loading Data with Minfi Method ]")
    message("----------------------------------")
    message("Loading data from ", directory)
    myDir <- directory
    # suppressWarnings(targets <- read.metharray.sheet(myDir)) # changed
    # new line
    suppressWarnings(targets <- read.metharray.sheet(myDir, pattern = sampleSheet.csv))
    rgSet <- read.metharray.exp(targets = targets, extended = TRUE,
                                force = force)
    if (arraytype == "EPIC")
      rgSet@annotation <- c(array = "IlluminaHumanMethylationEPIC",
                            annotation = "ilm10b4.hg19")
    sampleNames(rgSet) = rgSet[[1]]
    pd <- pData(rgSet)

    # new lines (if else if statement only)
    if (preproc == "Raw"){
      message("No preprocessing was done! This is the default setting for ChAMP")
      mset <- preprocessRaw(rgSet) # default line

    } else if (preproc == "Noob"){
      message("Normexp Convolution model is used for background correction! \n")
      mset <- preprocessNoob(rgSet, dyeMethod = dyeMethod, verbose = TRUE)
    }
    detP <- detectionP(rgSet)
    message("<< Read DataSet Success. >>\n")
    if (methValue == "B")
      tmp = getBeta(mset, "Illumina")
    else tmp = getM(mset)
    tmp[detP >= detPcut] <- NA
    message("The fraction of failed positions per sample\n \n            (You may need to delete samples with high proportion of failed probes\n): ")
    numfail <- matrix(colMeans(is.na(tmp)))
    rownames(numfail) <- colnames(detP)
    colnames(numfail) <- "Failed CpG Fraction."
    print(numfail)
    RemainSample <- which(numfail < SampleCutoff)
    if (any(numfail >= SampleCutoff))
      message("The detSamplecut parameter is : ", SampleCutoff,
              "\nSamples : ", paste(rownames(numfail)[which(numfail >=
                                                              SampleCutoff)], collapse = ","), " will be deleted.\n",
              "There are ", length(RemainSample), " samples left for analysis.\n")
    rgSet <- rgSet[, RemainSample]
    detP <- detP[, RemainSample]
    mset <- mset[, RemainSample]
    pd <- pd[RemainSample, ]
    tmp <- tmp[, RemainSample]
    if (filterDetP) {
      mset.f = mset[rowSums(is.na(tmp)) <= ProbeCutoff *
                      ncol(detP), ]
      if (ProbeCutoff == 0) {
        message("Filtering probes with a detection p-value above ",
                detPcut, " in one or more samples has removed ",
                dim(mset)[1] - dim(mset.f)[1], " probes from the analysis. If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples.")
      }
      else {
        message("Filtering probes with a detection p-value above ",
                detPcut, " in at least ", ProbeCutoff * 100,
                "% of samples has removed ", dim(mset)[1] -
                  dim(mset.f)[1], " probes from the analysis. If a large number of probes have been removed, ChAMP suggests you look at the failedSample file to identify potentially bad samples.")
      }
      mset = mset.f
      tmp <- tmp[rowSums(is.na(tmp)) <= ProbeCutoff *
                   ncol(detP), ]
      message("<< Filter DetP Done. >>\n")
    }
    if (sum(is.na(tmp)) == 0) {
      message("\nThere is no NA values in your matrix, there is no need to imputation.\n")
    }
    else {
      message("\nThere are ", sum(is.na(tmp)), " NA remain in filtered Data Set. Impute can be done for remain NAs, but not suitable for small number samples. For small Data Set (like only 20 samples), we suggest you set parameter ProbeCutoff as 0 in champ.load() here, which would remove all NA involved probe no matter how many samples of those probes are NA.\n")
    }
    if (autoimpute & sum(is.na(tmp)) > 0) {
      message("Impute will be conducted here for remain ",
              sum(is.na(tmp)), "  NAs. Note that if you don't do this, NA values would be kept in your data set. You may use champ.impute() function to do more complex imputation as well.")
      message("\nImpute function is working now, it may need couple minutes...")
      zz <- file("ImputeMessage.Rout", open = "wt")
      sink(zz)
      sink(zz, type = "message")
      tmp <- impute.knn(tmp, k = 5)$data
      sink(type = "message")
      sink()
      message("<< Imputation Done. >>\n")
    }
    if (filterBeads) {
      bc = mybeadcount(rgSet)
      bc2 = bc[rowSums(is.na(bc)) < beadCutoff * (ncol(bc)),
      ]
      mset.f2 = mset[featureNames(mset) %in% row.names(bc2),
      ]
      tmp <- tmp[rownames(tmp) %in% row.names(bc2), ]
      message("Filtering probes with a beadcount <3 in at least ",
              beadCutoff * 100, "% of samples, has removed ",
              dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
      mset = mset.f2
      message("<< Filter Beads Done. >>\n")
    }
    if (filterNoCG) {
      mset.f2 = dropMethylationLoci(mset, dropCH = T)
      tmp <- tmp[rownames(tmp) %in% featureNames(mset.f2),
      ]
      message("Filtering non-cg probes, has removed ",
              dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
      mset <- mset.f2
      message("<< Filter NoCG Done. >>\n")
    }
    if (filterSNPs) {
      if (arraytype == "450K") {
        if (is.null(population)) {
          message("Using general 450K SNP list for filtering.")
          data(hm450.manifest.hg19)
          maskname <- rownames(hm450.manifest.hg19)[which(hm450.manifest.hg19$MASK_general ==
                                                            TRUE)]
        }
        else if (!population %in% c("AFR", "EAS", "EUR",
                                    "SAS", "AMR", "GWD", "YRI", "TSI", "IBS",
                                    "CHS", "PUR", "JPT", "GIH", "CHB", "STU",
                                    "ITU", "LWK", "KHV", "FIN", "ESN", "CEU",
                                    "PJL", "ACB", "CLM", "CDX", "GBR", "BEB",
                                    "PEL", "MSL", "MXL", "ASW")) {
          message("Using general 450K SNP list for filtering.")
          data(hm450.manifest.hg19)
          maskname <- rownames(hm450.manifest.hg19)[which(hm450.manifest.hg19$MASK_general ==
                                                            TRUE)]
        }
        else {
          message("Using ", population, " specific 450K SNP list for filtering.")
          data(hm450.manifest.pop.hg19)
          maskname <- rownames(hm450.manifest.pop.hg19)[which(hm450.manifest.pop.hg19[,
                                                                                      paste("MASK_general_", population, sep = "")] ==
                                                                TRUE)]
        }
      }
      else {
        if (is.null(population)) {
          message("Using general EPIC SNP list for filtering.")
          data(EPIC.manifest.hg19)
          maskname <- rownames(EPIC.manifest.hg19)[which(EPIC.manifest.hg19$MASK_general ==
                                                           TRUE)]
        }
        else if (!population %in% c("AFR", "EAS", "EUR",
                                    "SAS", "AMR", "GWD", "YRI", "TSI", "IBS",
                                    "CHS", "PUR", "JPT", "GIH", "CHB", "STU",
                                    "ITU", "LWK", "KHV", "FIN", "ESN", "CEU",
                                    "PJL", "ACB", "CLM", "CDX", "GBR", "BEB",
                                    "PEL", "MSL", "MXL", "ASW")) {
          message("Using general EPIC SNP list for filtering.")
          data(EPIC.manifest.hg19)
          maskname <- rownames(EPIC.manifest.hg19)[which(EPIC.manifest.hg19$MASK_general ==
                                                           TRUE)]
        }
        else {
          message("Using ", population, " specific EPIC SNP list for filtering.")
          data(EPIC.manifest.pop.hg19)
          maskname <- rownames(EPIC.manifest.pop.hg19)[which(EPIC.manifest.pop.hg19[,
                                                                                    paste("MASK_general_", population, sep = "")] ==
                                                               TRUE)]
        }
      }
      mset.f2 = mset[!featureNames(mset) %in% maskname,
      ]
      tmp <- tmp[!rownames(tmp) %in% maskname, ]
      message("Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper, 2016, has removed ",
              dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
      mset = mset.f2
      message("<< Filter SNP Done. >>\n")
    }
    if (filterMultiHit) {
      data(multi.hit)
      mset.f2 = mset[!featureNames(mset) %in% multi.hit$TargetID,
      ]
      tmp <- tmp[!rownames(tmp) %in% multi.hit$TargetID,
      ]
      message("Filtering probes that align to multiple locations as identified in Nordlund et al, has removed ",
              dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
      mset = mset.f2
      message("<< Filter MultiHit Done. >>\n")
    }
    if (filterXY) {
      if (arraytype == "EPIC")
        data(probe.features.epic)
      else data(probe.features)
      autosomes = probe.features[!probe.features$CHR %in%
                                   c("X", "Y"), ]
      mset.f2 = mset[featureNames(mset) %in% row.names(autosomes),
      ]
      tmp <- tmp[rownames(tmp) %in% row.names(autosomes),
      ]
      message("Filtering probes on the X or Y chromosome has removed ",
              dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")
      mset = mset.f2
      message("<< Filter XY chromosome Done. >>\n")
    }
    message(paste(if (methValue == "B")
      "[Beta"
      else "[M", "value is selected as output.]\n"))
    beta.raw <- tmp
    intensity <- minfi::getMeth(mset) + minfi::getUnmeth(mset)
    detP <- detP[which(row.names(detP) %in% row.names(beta.raw)),
    ]
    if (min(beta.raw, na.rm = TRUE) <= 0)
      beta.raw[beta.raw <= 0] <- min(beta.raw[beta.raw >
                                                0])
    message("Zeros in your dataset have been replaced with smallest positive value.\n")
    if (max(beta.raw, na.rm = TRUE) >= 0)
      beta.raw[beta.raw >= 1] <- max(beta.raw[beta.raw <
                                                1])
    message("One in your dataset have been replaced with largest value below 1.\n")
    message("The analysis will be proceed with ", dim(beta.raw)[1],
            " probes and ", dim(beta.raw)[2], " samples.\n")
    message("Current Data Set contains ", sum(is.na(beta.raw)),
            " NA in ", if (methValue == "B")
              "[Beta]"
            else "[M]", " Matrix.\n")
    message("[<<<<< ChAMP.LOAD END >>>>>>]")
    message("[===========================]")
    message("[You may want to process champ.QC() next.]\n")
    return(list(mset = mset, rgSet = rgSet, pd = pd, intensity = intensity,
                beta = beta.raw, detP = detP))
  }
  else {
    message("\n[ Loading Data with ChAMP Method ]")
    message("----------------------------------")
    message("Note that ChAMP method will NOT return rgSet or mset, they object defined by minfi. Which means, if you use ChAMP method to load data, you can not use SWAN or FunctionNormliazation method in champ.norm() (you can use BMIQ or PBC still). But All other function should not be influenced.\n")
    myImport <- champ.import(directory, arraytype = arraytype)
    if (methValue == "B")
      myLoad <- champ.filter(beta = myImport$beta, M = NULL,
                             pd = myImport$pd, intensity = myImport$intensity,
                             Meth = NULL, UnMeth = NULL, detP = myImport$detP,
                             beadcount = myImport$beadcount, autoimpute = autoimpute,
                             filterDetP = filterDetP, ProbeCutoff = ProbeCutoff,
                             SampleCutoff = SampleCutoff, detPcut = detPcut,
                             filterBeads = filterBeads, beadCutoff = beadCutoff,
                             filterNoCG = filterNoCG, filterSNPs = filterSNPs,
                             population = population, filterMultiHit = filterMultiHit,
                             filterXY = filterXY, arraytype = arraytype)
    else myLoad <- champ.filter(beta = NULL, M = myImport$M,
                                pd = myImport$pd, intensity = myImport$intensity,
                                Meth = NULL, UnMeth = NULL, detP = myImport$detP,
                                beadcount = myImport$beadcount, autoimpute = autoimpute,
                                filterDetP = filterDetP, ProbeCutoff = ProbeCutoff,
                                SampleCutoff = SampleCutoff, detPcut = detPcut,
                                filterBeads = filterBeads, beadCutoff = beadCutoff,
                                filterNoCG = filterNoCG, filterSNPs = filterSNPs,
                                population = population, filterMultiHit = filterMultiHit,
                                filterXY = filterXY, arraytype = arraytype)
    message("[<<<<< ChAMP.LOAD END >>>>>>]")
    message("[===========================]")
    message("[You may want to process champ.QC() next.]\n")
    return(myLoad)
  }
}
