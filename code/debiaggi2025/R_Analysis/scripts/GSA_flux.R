# Libraries (merging data)
library(reshape2)
library(readxl)
library(readr)

# Libraries (from SILAC to absolute quantification)
library(ggplot2)
library(biomaRt)
library(knitr)
library(plyr);
library(dplyr)
library(ggrepel)

# Libraries for using Piano package
library(piano)
library(biomartr)
library(limma)
library(gtools)

rm(list = ls())#remove everything from environment

GSA_list = list("ecCARe_vs_ecPARe",
                "ecOBEe_vs_ecPARe",
                "ecCARn_vs_ecPARn",
                "ecOBEn_vs_ecPARn")

# Load GO terms
GO <- read.csv("data/processed/rxnGO_split.csv", header = TRUE, sep = ",", quote = "\"")
myGsc <- loadGSC(GO)

for (strain_pair in GSA_list){
  
  print(paste("Working on: ", strain_pair, "fluxes"))
  
  # Gene Set Analysis ------------------------------------- doi:10.1093/nar/gkt111.
  # Load df and adj-pvalue data
  
  GSA_source_path <- file.path("data/processed",paste0("flux_analysis_results_",strain_pair,".csv"))
  
  MSdata <- read.csv(GSA_source_path,header=TRUE, sep = ";", quote = "\"")
  
  # Replace Inf and -Inf with NA
  MSdata[sapply(MSdata, is.infinite)] <- NA
  
  # Remove rows with NA values
  MSdata <- na.omit(MSdata)  # Removes rows with any NA in the dataframe
  
  MSdata_p <- as.matrix(MSdata[,3])
  rownames(MSdata_p) <- MSdata[,1]
  MSdata_l <- as.matrix(MSdata[,2])
  rownames(MSdata_l) <- MSdata[,1]
  
  gsaRes_MSdata <- runGSA(MSdata_p,MSdata_l,gsc=myGsc,gsSizeLim=c(3,300))
  
  GSA_output_indiv_path <- file.path("data/processed",paste0("gsaResTab_flux_",strain_pair,".csv"))
  GSA_output_indiv <- GSAsummaryTable(gsaRes_MSdata)
  write.csv(GSA_output_indiv, file = GSA_output_indiv_path, row.names = FALSE)
}