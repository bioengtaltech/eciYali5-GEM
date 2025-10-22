# Load necessary libraries
library(dplyr)
library(ggplot2)
library(PCAtools)
library(scales)
library(svglite)

# Load proteomics data
proteomics_path <- "data/processed/intensitiesTable.tsv"
proteomics_data <- read.delim(proteomics_path, header = TRUE, sep = "\t")

# Extract protein names and intensity values
protein_names <- proteomics_data[[1]]             # First column: protein names
intensity_values <- proteomics_data[, -1]         # Exclude protein names for PCA

# Remove rows with NaN values
complete_data <- na.omit(intensity_values)

# Create metadata from column names
metadata <- data.frame(Sample = colnames(complete_data))
rownames(metadata) <- metadata$Sample

# Extract Group_ID (e.g., PAR, CAR, OBE) and Condition (e/n) from names like "PARe.1"
metadata$Group_ID <- sub("^([A-Z]+)[en]\\..*", "\\1", metadata$Sample)
metadata$Condition <- sub("^[A-Z]+([en])\\..*", "\\1", metadata$Sample)

# Order Group_ID and rename Condition for clarity
metadata$Group_ID <- factor(metadata$Group_ID, levels = c("PAR", "CAR", "OBE"))
metadata$Condition <- factor(metadata$Condition, levels = c("e", "n"))
levels(metadata$Condition) <- c("Exp", "Nlim")

# Perform PCA
pca_result <- pca(complete_data, metadata = metadata, center = TRUE, scale = TRUE)

# Prepare data for plotting
plotdata <- data.frame(
  x = pca_result$rotated[, 1],
  y = pca_result$rotated[, 2],
  Group_ID = pca_result$metadata$Group_ID,
  Condition = pca_result$metadata$Condition
)

# Compute symmetric axis limits
max_axis_limit <- max(abs(plotdata$x), abs(plotdata$y))

# Define output file and save as SVG
output_file <- "results/figures/pca_biplot_fixed_axes.svg"
svglite(output_file, width = 6.75, height = 6.75)

# Create PCA biplot
p <- ggplot(plotdata, aes(x = x, y = y, color = Group_ID, shape = Condition)) +
  geom_point(size = 3, stroke = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  scale_shape_manual(values = c("Exp" = 19, "Nlim" = 1), name = "Growth Phase") +
  scale_color_discrete(name = "Strain") +
  coord_fixed(ratio = 1) +
  scale_x_continuous(
    limits = c(-max_axis_limit, max_axis_limit),
    name = paste0("PC1 (", round(pca_result$variance[1], 2), "%)")
  ) +
  scale_y_continuous(
    limits = c(-max_axis_limit, max_axis_limit),
    name = paste0("PC2 (", round(pca_result$variance[2], 2), "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_blank()
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2, override.aes = list(stroke = 1.5, size = 3))
  )

# Print and save
print(p)
ggsave(output_file, plot = p, width = 6.75, height = 6.75)
