# Load necessary library
library(dplyr)

# Load your .tsv files
gene_intensities_path <- "data/raw/ggMatrix.tsv"
ggmatrix <- read.delim(gene_intensities_path, header = TRUE, sep = "\t")

uniprot_IDs_YaliCLIB_path <- "data/raw/uniprotkb_taxonomy_id_284591.tsv"
YaliCLIB_uniprot_IDs <- read.delim(uniprot_IDs_YaliCLIB_path, header = TRUE, sep = "\t")

uniIDs_Yali_path <- "data/raw/uniprotkb_taxonomy_id_4952.tsv"
Yali_uniIDs <- read.delim(uniIDs_Yali_path, header = TRUE, sep = "\t")

# Split the 'Genes' column into separate gene names
ggmatrix <- ggmatrix %>%
  mutate(Genes_List = strsplit(as.character(Genes), ";"))

# Split Gene Names into a list of names for both UniProt datasets
Yali_uniIDs <- Yali_uniIDs %>%
  mutate(Gene_Names_List = strsplit(as.character(Gene.Names), " "))

YaliCLIB_uniprot_IDs <- YaliCLIB_uniprot_IDs %>%
  mutate(Gene_Names_List = strsplit(as.character(Gene.Names), " "))

# Define the function to get a single UniProt entry based on gene names
get_uniprot_entry <- function(genes_list, uniprot_df) {
  # Ensure genes_list is a vector of gene names
  genes_vector <- unlist(strsplit(genes_list, ", "))
  
  # Check for the first matching entry in the UniProt data frame
  matching_entry <- uniprot_df %>%
    filter(sapply(Gene_Names_List, function(names) any(genes_vector %in% names))) %>%
    pull(Entry) %>%
    .[1]  # Take only the first match
  
  if (!is.na(matching_entry)) {
    return(matching_entry)
  } else {
    return(NA)
  }
}

# Apply the lookup function first using YaliCLIB_uniprot_IDs, and then fill missing entries with Yali_uniIDs
ggmatrix <- ggmatrix %>%
  mutate(UniProt_Entry = sapply(Genes_List, function(x) {
    entry <- get_uniprot_entry(x, YaliCLIB_uniprot_IDs)
    if (is.na(entry)) {
      entry <- get_uniprot_entry(x, Yali_uniIDs)  # If missing, try the broader taxonomy group
    }
    return(entry)
  }))

# Order the intensity columns by biological prefix pattern
desired_order <- c("PARe\\.", "PARn\\.", "CARe\\.", "CARn\\.", "OBEe\\.", "OBEn\\.")

# Select and order intensity columns according to the defined prefix order
intensity_columns <- unlist(lapply(desired_order, function(pattern) {
  grep(paste0("^", pattern), colnames(ggmatrix), value = TRUE)
}))

# Create the new data frame with UniProt_Entry and ordered intensity columns
intensitiesTable <- ggmatrix %>%
  select(UniProt_Entry, all_of(intensity_columns))  # Keep only the UniProt_Entry and ordered intensity columns

# Remove proteins with all intensity values missing (unused strains)
intensitiesTable <- intensitiesTable %>%
  rowwise() %>%
  mutate(all_intensities_NA = all(is.na(c_across(all_of(intensity_columns))))) %>%
  ungroup() %>%
  filter(!all_intensities_NA) %>%
  select(-all_intensities_NA)

# Write the processed table to file
output_path <- "data/processed/intensitiesTable.tsv"
write.table(intensitiesTable, file = output_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

