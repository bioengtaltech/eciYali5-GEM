# Load necessary libraries
library(dplyr)
library(readr)

# Load data using load_ggMatrix.R (make sure it's in the working directory or provide a path)
source("scripts/load_ggMatrix.R")  # This generates `intensitiesTable` and loads `uniprotData`

# Load the calcTPA function from its script
source("scripts/calcTPA.R")  # Makes calcTPA available for use

# Load the trypsin_digestion function from its script
source("scripts/trypsin_digestion.R")  # Makes trypsin_digestion available for normalization

# Load proteinContents.tsv
proteinContents_path <- "data/raw/proteinContents.tsv"
proteinContents <- read.delim(proteinContents_path, header = TRUE, sep = "\t")

# Normalization factor (adjustable)
normFactor <- 0.99

# Function to calculate protein abundances for all samples in `intensitiesTable`
calculate_abundances <- function(intensitiesTable, uniprotData, proteinContents, normFactor) {
  # Run calcTPA to get protein abundances
  proteinAbundances <- calcTPA(intensitiesTable, uniprotData, normalize = TRUE)
  
  # Prepare an output data frame to store results for each sample
  abundancesForGECKO <- data.frame(UniProt_Entry = proteinAbundances$UniProt_Entry)
  
  # Identify columns representing sample intensities (by biological prefixes)
  intensity_prefixes <- c("PARe\\.", "PARn\\.", "CARe\\.", "CARn\\.", "OBEe\\.", "OBEn\\.")
  intensity_cols <- unlist(lapply(intensity_prefixes, function(p) {
    grep(paste0("^", p), colnames(intensitiesTable), value = TRUE)
  }))
  
  # Loop through each sample column to calculate abundances for GECKO
  for (sample_col in intensity_cols) {
    # Extract strain name (e.g., PARe, PARn, etc.)
    strain_name <- sub("\\..*", "", sample_col)
    
    # Match Ptot from proteinContents.tsv based on strain name
    Ptot <- proteinContents$Ptot[match(strain_name, proteinContents$Strain)]
    
    # Check if Ptot is found for the strain
    if (is.na(Ptot)) {
      warning(paste("Ptot not found for strain:", strain_name))
      next
    }
    
    # Calculate abundForGECKO: PredictedAbundance * 10^3 * Ptot * normFactor
    abundancesForGECKO[[sample_col]] <- 
      proteinAbundances[[sample_col]] * 10^3 * Ptot * normFactor
  }
  
  # Return both data frames
  return(list(abundancesForGECKO = abundancesForGECKO,
              proteinAbundances = proteinAbundances))
}

# Run the function and store results
result_list <- calculate_abundances(intensitiesTable, Yali_uniIDs, proteinContents, normFactor)

# Extract `abundForGECKO` and `proteinAbundances` from the returned list
abundForGECKO <- result_list$abundancesForGECKO
proteinAbundances <- result_list$proteinAbundances

# Save `abundForGECKO` as a .tsv file
output_path <- "data/processed/proteomics4GECKO.tsv"
write.table(abundForGECKO, file = output_path, sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# Save `proteinAbundances` as a separate .tsv file
proteinAbundances_path <- "data/processed/proteinAbundances.tsv"
write.table(proteinAbundances, file = proteinAbundances_path, sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)