## Script automatically generated on Wed Feb 17 09:44:46 2021

library(patRoon)
library(xcms)

# -------------------------
# initialization
# -------------------------

options(patRoon.path.MetFragCL = "D:/MACHINE LEARNING/projects/xcms/MetFrag2.4.5-CL.jar")

workPath <- "D:/MACHINE LEARNING/projects/xcms"
setwd(workPath)

anaInfo <- generateAnalysisInfo(paths = c("D:/MACHINE LEARNING/projects/xcms/asd"),
                                groups = c("190524_Int_002", "190524_Int_012", "190524_Int_019"),
                                blanks = c("", "", ""))

# Set to FALSE to skip data pre-treatment
doDataPretreatment <- TRUE
if (doDataPretreatment)
{
    recalibrarateDAFiles(anaInfo) 
} 


# -------------------------
# features
# -------------------------


# Find all features.
# NOTE: see XCMS manual for many more options
fList <- findFeatures(anaInfo, "xcms3", param = xcms::CentWaveParam(peakwidth = c(5, 15)))

# Group and align features between analysis
fGroups <- groupFeatures(fList, "xcms3", rtalign = TRUE,
                         groupParam = xcms::PeakDensityParam(sampleGroups = analysisInfo(fList)$group),
                         retAlignParam = xcms::ObiwarpParam())

# Basic rule based filtering
fGroups <- filter(fGroups, preAbsMinIntensity = 100, absMinIntensity = 10000,
                  relMinReplicateAbundance = 1, maxReplicateIntRSD = 0.75,
                  blankThreshold = 5, removeBlanks = TRUE,
                  retentionRange = NULL, mzRange = NULL)


# -------------------------
# annotation
# -------------------------


# Retrieve MS peak lists
avgPListParams <- getDefAvgPListParams(clusterMzWindow = 0.005)
mslists <- generateMSPeakLists(fGroups, "mzr", maxMSRtWindow = NULL, precursorMzWindow = NULL,
                               avgFeatParams = avgPListParams, avgFGroupParams = avgPListParams)
# uncomment and configure for extra filtering of MS peak lists
mslists <- filter(mslists, absMSIntThr = NULL, absMSMSIntThr = NULL, relMSIntThr = NULL,
                 relMSMSIntThr = NULL, topMSPeaks = NULL, topMSMSPeaks = NULL,
                 deIsotopeMS = FALSE, deIsotopeMSMS = FALSE)

# Calculate formula candidates
formulas <- generateFormulas(fGroups, "genform", mslists, relMzDev = 5,
                             adduct = "[M+H]+", elements = "CHNOP",
                             calculateFeatures = TRUE, featThreshold = 0.75)

# Find compound structure candidates
compounds <- generateCompounds(fGroups, mslists, "metfrag", method = "R", dbRelMzDev = 5,
                               fragRelMzDev = 5, fragAbsMzDev = 0.002,
                               adduct = "[M+H]+", database = "pubchem", maxCandidatesToStop = 2500)

compounds <- addFormulaScoring(compounds, formulas, TRUE) 


# Perform automatic generation of components
components <- generateComponents(fGroups, "ramclustr", ionization = "positive")


# -------------------------
# reporting
# -------------------------


reportCSV(fGroups, path = "report", reportFeatures = FALSE, formulas = formulas,
          compounds = compounds, compoundsNormalizeScores = "max",
          components = components)

reportHTML(fGroups, path = "report", reportPlots = c("chord", "venn", "upset", "eics", "formulas"),
           formulas = formulas, compounds = compounds, compoundsNormalizeScores = "max",
           components = components, MSPeakLists = mslists,
           selfContained = FALSE, openReport = TRUE)

