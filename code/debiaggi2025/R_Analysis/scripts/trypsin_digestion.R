trypsin_digestion <- function(sequence, minLength, maxLength) {
  # Load the stringr package
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("The stringr package is required but is not installed. Please install it using install.packages('stringr').")
  }
  
  # Validate the input
  if (!is.character(sequence)) {
    stop("The input sequence must be a string.")
  }
  
  # Cleavage site for trypsin (C-terminal side of arginine or lysine, except when followed by proline)
  cleavageSite <- "(?<=R|K)(?!P)"
  
  # Find the cleavage sites using stringr's str_locate_all
  cleavageIndices <- stringr::str_locate_all(sequence, cleavageSite)[[1]][, "end"]
  
  # Preallocate the list for peptides
  peptides <- list()
  
  # Extract the peptides
  startIndex <- 1
  for (i in seq_along(cleavageIndices)) {
    endIndex <- cleavageIndices[i]
    peptides[[i]] <- substr(sequence, startIndex, endIndex)
    startIndex <- endIndex + 1
  }
  
  # Add the last peptide (remaining sequence)
  peptides[[length(cleavageIndices) + 1]] <- substr(sequence, startIndex, nchar(sequence))
  
  # Filter peptides by length
  valid_peptides <- sapply(peptides, function(pep) {
    nchar(pep) >= minLength && nchar(pep) <= maxLength
  })
  
  peptides <- peptides[valid_peptides]
  
  # Count the number of predicted peptides
  numPeptides <- length(peptides)
  
  return(list(peptides = peptides, numPeptides = numPeptides))
}