# Load necessary libraries
library(dplyr)
library(tidyr)
library(stringr)

# Load the data
geneGOdata <- read.delim("data/raw/yali_genes_GOterms.tsv", header = TRUE, sep = "\t")
rxnGenedata <- read.csv("data/raw/rxnGoTermsTable.csv", header = TRUE, sep = ";", quote = "\"")

# Prepare geneGOdata for matching
geneGOdata <- geneGOdata %>%
  mutate(Gene.Names = strsplit(as.character(Gene.Names), " ")) %>% # Split Gene.Names into lists
  unnest_longer(Gene.Names) # Expand lists into separate rows

# Process rxnGenedata
rxnGenedata <- rxnGenedata %>%
  mutate(genes = gsub("\\(|\\)", "", as.character(genes)), # Remove parentheses
         genes = strsplit(genes, " (and|or) ")) %>% # Split genes on 'and' or 'or'
  unnest_longer(genes) %>% # Expand lists into separate rows
  mutate(genes = trimws(genes)) # Trim whitespace from gene names

# Join data and extract GO terms
# Perform the join and handle missing matches
rxnGO <- rxnGenedata %>%
  left_join(geneGOdata, by = c("genes" = "Gene.Names")) %>% # Match genes to Gene.Names
  filter(!is.na(Gene.Ontology..biological.process.)) %>% # Exclude unmatched rows
  group_by(rxns) %>% # Group by Entry
  summarise(GO_terms = paste(unique(Gene.Ontology..biological.process.), collapse = "; ")) # Aggregate unique GO terms

# Split rxnGO so each rxn-GO_term pair has its own row
rxnGO_split <- rxnGO %>%
  separate_rows(GO_terms, sep = "; ") %>% # Split GO_terms into separate rows
  mutate(GO_terms = trimws(GO_terms)) # Trim any extra whitespace

# Remove rows with blank or missing GO_terms
rxnGO_split <- rxnGO_split %>%
  filter(GO_terms != "", !is.na(GO_terms)) # Remove blank or NA GO_terms

# Remove duplicated rows based on both rxn and GO_terms
rxnGO_split <- rxnGO_split %>%
  distinct(rxns, GO_terms, .keep_all = TRUE) # Remove duplicates

# View the cleaned result
write.csv(rxnGO_split, "data/processed/rxnGO_split.csv", row.names = FALSE)