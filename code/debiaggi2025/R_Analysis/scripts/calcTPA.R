# calcTPA function
calcTPA <- function(proteinIntensities, uniprotData, normalize = FALSE) {
  
  # Match protein intensities with UniProt data
  matchingIndices <- match(proteinIntensities$UniProt_Entry, uniprotData$Entry)
  proteinIntensities$MW <- uniprotData$Mass[matchingIndices]
  
  # Prepare output: UniProt entries in the first column
  result <- data.frame(UniProt_Entry = proteinIntensities$UniProt_Entry)
  
  # Identify columns representing intensity data (by biological prefixes)
  intensity_prefixes <- c("PARe\\.", "PARn\\.", "CARe\\.", "CARn\\.", "OBEe\\.", "OBEn\\.")
  intensity_cols <- unlist(lapply(intensity_prefixes, function(p) {
    grep(paste0("^", p), colnames(proteinIntensities), value = TRUE)
  }))
  
  # Loop through each intensity column (sample) to calculate protein abundances
  for (sample_col in intensity_cols) {
    sample_intensities <- proteinIntensities[[sample_col]]
    
    if (any(is.na(proteinIntensities$MW))) {
      stop("Some molecular weights (MW) are missing. Check the UniProt data and protein intensities for mismatches.")
    }
    
    # Non-normalized case
    if (!normalize) {
      totalIntensity <- sum(proteinIntensities$MW * sample_intensities, na.rm = TRUE)
      predictedAbundance <- sample_intensities / totalIntensity
      predictedAbundance <- predictedAbundance * proteinIntensities$MW
    } else {
      # Normalized case with theoretical peptide correction
      proteinIntensities$theorPep <- sapply(uniprotData$Sequence[matchingIndices], function(seq) {
        digestion <- trypsin_digestion(seq, 6, 30)
        digestion$numPeptides
      })
      
      totalIntensity <- sum(proteinIntensities$MW * sample_intensities / proteinIntensities$theorPep, na.rm = TRUE)
      predictedAbundance <- (sample_intensities / proteinIntensities$theorPep) / totalIntensity
      predictedAbundance <- predictedAbundance * proteinIntensities$MW
    }
    
    # Add the calculated abundances for this sample as a new column in the result
    result[[sample_col]] <- predictedAbundance
  }
  
  # Return the result containing UniProt entries and calculated abundances for all samples
  return(result)
}