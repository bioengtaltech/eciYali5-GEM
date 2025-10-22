# Load necessary libraries
library(dplyr)

# Load the data
data <- read.delim("data/raw/QuickGO-annotations-Yali.tsv", header = TRUE, sep = "\t")

# Transform to the desired format
transformed_data <- data %>%
  mutate(GOgroups = paste(GO.NAME, " [", GO.TERM, "]", sep = "")) %>%  # Combine GO.NAME and GO.TERM
  select(Entry = GENE.PRODUCT.ID, time = GO.ASPECT, GOgroups)         # Rename and select relevant columns

# Save transformed_data to a .csv file
write.csv(transformed_data, "data/processed/QuickGO-annotations-Yali.csv", row.names = FALSE)

# Confirm the file was saved
cat("Transformed data saved to 'data/processed/QuickGO-annotations-Yali.csv'")