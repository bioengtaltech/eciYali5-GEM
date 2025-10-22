# Load required packages
library(dplyr)
library(tibble)
library(ggplot2)
library(EnhancedVolcano)

# Load your .tsv file
protein_abundances_path <- "data/processed/proteinAbundances.tsv"
proteinAbundances <- read.delim(protein_abundances_path, header = TRUE, sep = "\t")

# Select columns of interest for the two samples

comparison_list <- list(
  list(group1_name = "PARe", group2_name = "CARe", comparison_label = "CARe_vs_PARe",
       exclude_group1 = NULL, exclude_group2 = NULL),
  list(group1_name = "PARe", group2_name = "OBEe", comparison_label = "OBEe_vs_PARe",
       exclude_group1 = NULL, exclude_group2 = NULL),
  list(group1_name = "PARn", group2_name = "CARn", comparison_label = "CARn_vs_PARn",
       exclude_group1 = NULL, exclude_group2 = c("CARn.3")),
  list(group1_name = "PARn", group2_name = "OBEn", comparison_label = "OBEn_vs_PARn",
       exclude_group1 = NULL, exclude_group2 = NULL)
)

all_comparison_data_list <- list()

for (comparison in comparison_list) {
  
  print(paste("Processing comparison:", comparison$comparison_label))
  
  # sample column determination
  group1_name <- comparison$group1_name
  group2_name <- comparison$group2_name
  
  sample1_cols <- grep(paste0("^", group1_name,"\\."), colnames(proteinAbundances), value = TRUE)
  sample2_cols <- grep(paste0("^", group2_name,"\\."), colnames(proteinAbundances), value = TRUE)
  
  sample1_cols <- sample1_cols[!(sample1_cols %in% comparison$exclude_group1)]
  sample2_cols <- sample2_cols[!(sample2_cols %in% comparison$exclude_group2)]
  
  # Remove rows with any missing values in the selected columns using `if_all`
  proteinAbundances_clean <- proteinAbundances %>%
    filter(!is.na(UniProt_Entry)) %>%
    filter(if_all(all_of(c(sample1_cols, sample2_cols)), ~ !is.na(.)))
  
  # Calculate sample means and log2 fold changes
  comparison_df <- proteinAbundances_clean %>%
    rowwise() %>%
    mutate(
      sample1_mean = mean(c_across(all_of(sample1_cols)), na.rm = TRUE),
      sample2_mean = mean(c_across(all_of(sample2_cols)), na.rm = TRUE),
      log2FC = log2(sample2_mean / sample1_mean)
    ) %>%
    ungroup()  # Important: Ungroup the dataframe to avoid errors when adding vector columns
  
  # Conduct a Student's T-test for each protein and obtain p-values
  p_values <- apply(proteinAbundances_clean[, c(sample1_cols, sample2_cols)], 1, function(x) {
    sample1 <- x[1:length(sample1_cols)]
    sample2 <- x[(length(sample1_cols) + 1):(length(x))]
    
    if (sum(!is.na(sample1)) >= 2 && sum(!is.na(sample2)) >= 2) {
      return(t.test(sample1, sample2, na.action = na.omit)$p.value)
    } else {
      return(NA)  # Return NA if there are not enough observations in either group
    }
  })
  
  # Create a new data frame by combining comparison_df with p-values
  comparison_df <- comparison_df %>%
    mutate(p_value = p_values)  # Add p-values as a column
  
  # Apply multiple hypothesis correction (Benjamini-Hochberg)
  comparison_df <- comparison_df %>%
    mutate(adj_p_value = p.adjust(p_value, method = "BH"))  # Adjusted p-values using BH correction
  
  # Create a new data frame with only the columns you want: UniProt_Entry, log2FC, and adj_p_value
  final_comparison_df <- comparison_df %>%
    select(UniProt_Entry, log2FC, adj_p_value)
  
  log2FC_col_name <- paste0("log2FC_", comparison$comparison_label)
  adj_p_value_col_name <- paste0("adj_p_value_", comparison$comparison_label)
  
  colnames(final_comparison_df)[colnames(final_comparison_df) == "log2FC"] <- log2FC_col_name
  colnames(final_comparison_df)[colnames(final_comparison_df) == "adj_p_value"] <- adj_p_value_col_name
  
  all_comparison_data_list <- c(all_comparison_data_list, list(final_comparison_df))
  
  ## Uncomment to print individual pairwise DE analysis
  
  output_csv_path <- file.path("data/processed", paste0("DE_analysis_", comparison$comparison_label, ".csv"))
  write.csv(final_comparison_df, file = output_csv_path, row.names = FALSE)
  print(paste("Saved comparison data to:", output_csv_path))
  
  ## EXPLORATORY - Volcano plots (not used in publication)
  # Generate the volcano plot using EnhancedVolcano without labels
  
  # Calculate limits for x-axis and y-axis based on the data
  x_min <- min(comparison_df$log2FC, na.rm = TRUE) - 0.5   # Adding a small buffer
  x_max <- max(comparison_df$log2FC, na.rm = TRUE) + 0.5
  y_max <- max(-log10(comparison_df$adj_p_value), na.rm = TRUE) + 0.5  # Buffer for y-axis
  
  # Generate the volcano plot with adjusted x and y axis limits
  volcano_plot <- EnhancedVolcano(
    toptable = comparison_df,              
    lab = "",                               # No labels for data points,
    x = 'log2FC',                           # Column name for log2 fold change
    y = 'adj_p_value',                      # Column name for adjusted p-values
    xlab = bquote(~log[2]~"fold change"),   # X-axis label
    ylab = bquote(~-log[10]~P[adjusted]),   # Y-axis label
    title = NULL,                           # No title
    subtitle = NULL,                        # No subtitle
    pCutoff = 0.05,                         # Significance cut-off for adjusted p-values
    pCutoffCol = 'adj_p_value',             # Column used for p-value cut-off
    FCcutoff = 1.5,                         # Log2 fold-change cut-off
    cutoffLineWidth = 0.4,                  # Width of the cut-off lines
    pointSize = 3.0,                        # Size of the points representing proteins
    colAlpha = 0.5,                         # Transparency for plot points
    legendPosition = 'top',                 # Position of the legend
    legendLabels = c("NS", "log2 FC", "p-value", "p-value & log2 FC"), # Custom legend labels
    legendLabSize = 14,                     # Size of legend labels
    axisLabSize = 14,
    legendIconSize = 4.0,                   # Size of legend icons
    drawConnectors = FALSE,                 # Disable connectors since labels are off
    gridlines.major = FALSE,                # Remove major gridlines
    gridlines.minor = FALSE,                # Remove minor gridlines
    border = 'full',                        # Full border around the plot
    borderWidth = 0.5,                      # Width of the border
    xlim = c(x_min, x_max),                 # Custom x-axis limits
    ylim = c(0, y_max),                     # Custom y-axis limits
    caption = ""                            # No caption
  ) +
    # Adjust the axis line thickness to 0.5
    theme(
      axis.line = element_line(linewidth = 0.5),  # Set both x and y axis line thickness to 0.5
      axis.ticks = element_line(linewidth = .5)  # Set tick mark color and thickness
    )
  
  # Display the plot in the RStudio Plots pane
  print(volcano_plot)
  
  # Export the plot to an .svg file with specified dimensions
  output_svg_path <- file.path("results/figures", paste0("volcano_", comparison$comparison_label, ".svg"))
  ggsave(output_svg_path, plot = volcano_plot, device = "svg", width = 6.75, height = 6.75, units = "in")
  print(paste("Saved volcano plot to:", output_svg_path))
  
}

final_combined_table = Reduce(function(x, y) {full_join(x, y, by = "UniProt_Entry")}, all_comparison_data_list)

combinedTable_path <- file.path("data/processed", paste0("Table_S1_DE_analysis_all_pairs", ".csv"))
write.csv(final_combined_table, file = combinedTable_path, row.names = FALSE)
print(paste("Saved comparison data to:", combinedTable_path))