# Libraries (merging data)
library(reshape2)
library(readxl)
library(readr)

# Libraries (from SILAC to absolute quantification)
library(ggplot2)
library(biomaRt)
library(knitr)
library(plyr)
library(dplyr)
library(ggrepel)

# Libraries for using Piano package
library(piano)
library(biomartr)
library(limma)
library(gtools)

rm(list = ls()) # remove everything from environment

GSA_list = list(
  "CARe_vs_PARe",
  "OBEe_vs_PARe",
  "CARn_vs_PARn",
  "OBEn_vs_PARn",
  "PARn_vs_PARe",
  "CARn_vs_CARe",
  "OBEn_vs_OBEe"
)

# Load GO terms
GOraw <- read.csv("data/processed/QuickGO-annotations-Yali.csv", header = TRUE, sep = ",", quote = "\"")
GO <- GOraw[, c("Entry", "GOgroups")]

myGsc <- loadGSC(GO)

for (strain_pair in GSA_list) {
  
  print(paste("Working on:", strain_pair))
  
  GSA_source_path <- file.path("data/processed", paste0("DE_analysis_", strain_pair, ".csv"))
  
  # Gene Set Analysis ------------------------------------- doi:10.1093/nar/gkt111.
  # Load df and adj-pvalue data
  MSdata <- read.csv(GSA_source_path, header = TRUE)
  
  MSdata_p <- as.matrix(MSdata[, 3])
  rownames(MSdata_p) <- MSdata[, 1]
  MSdata_l <- as.matrix(MSdata[, 2])
  rownames(MSdata_l) <- MSdata[, 1]
  
  gsaRes_MSdata <- runGSA(MSdata_p, MSdata_l, gsc = myGsc, gsSizeLim = c(5, 300))
  
  GSA_output_indiv_path <- file.path("data/processed", paste0("gsaResTab_", strain_pair, ".csv"))
  GSA_output_indiv <- GSAsummaryTable(gsaRes_MSdata)
  write.csv(GSA_output_indiv, file = GSA_output_indiv_path, row.names = FALSE)
}