# ============================================================
# Libraries
# ============================================================
library(tidyverse)
library(ggrepel)
library(svglite)

# ============================================================
# Read Protein Abundances Table 
# ============================================================
intensities_df <- read_tsv("data/processed/intensitiesTable.tsv")

# ============================================================
# Define Groups and Samples to Exclude 
# ============================================================
group_definitions <- list(
  list(group_name = "PARe", exclude = NULL),
  list(group_name = "PARn", exclude = NULL),
  list(group_name = "CARe", exclude = NULL),
  list(group_name = "CARn", exclude = c("CARn.3")),
  list(group_name = "OBEe", exclude = NULL),
  list(group_name = "OBEn", exclude = NULL)
)

# ============================================================
# Compute Mean Intensities per Group 
# ============================================================
averaged_intensities <- intensities_df %>% select(UniProt_Entry)
for (group in group_definitions) {
  group_name <- group$group_name
  exclude_cols <- group$exclude
  group_cols <- grep(paste0("^", group_name, "\\."), colnames(intensities_df), value = TRUE)
  cols_to_average <- setdiff(group_cols, exclude_cols)
  mean_col_name <- paste0(group_name, "_mean")
  averaged_intensities <- averaged_intensities %>%
    mutate(!!mean_col_name := rowMeans(intensities_df[cols_to_average], na.rm = TRUE))
}

# ============================================================
# Highlighted Protein Sets
# ============================================================
set_PAR <- c("Q6C8T9", "Q6C704", "Q6C9V5", "Q9UVF4")
set_CAR <- c("Q6C8T9", "Q6C704", "A0A168PH23", "Q9UUQ6")
set_OBE <- c("Q6C9V5", "Q9UVF4")
highlight_definitions <- list(
  PARe = set_PAR, PARn = set_PAR,
  CARe = set_CAR, CARn = set_CAR,
  OBEe = set_OBE, OBEn = set_OBE
)

# ============================================================
# Protein Name Dictionary
# ============================================================
protein_names <- c(
  "Q6C8T9" = "GGS1", "Q6C704" = "HMGR", "Q6C9V5" = "DGA2",
  "Q9UVF4" = "GPD1", "Q6C0V2" = "Enolase", "A0A168PH23" = "CarB",
  "Q9UUQ6" = "CarRP"
)

# ============================================================
# Color Palette
# ============================================================
protein_colors <- c(
  "GGS1"    = "#0072B2", "HMGR"    = "#D55E00", "DGA2" = "#009E73",
  "GPD1"    = "#CC79A7", "Enolase" = "#56B4E9", "CarB" = "#E69F00",
  "CarRP"   = "#999999"
)

# ============================================================
# Custom Theme
# ============================================================
theme_pub <- function(base_size = 20, base_family = "Helvetica") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      text = element_text(size = base_size),
      axis.title.x = element_text(size = base_size),
      axis.title.y = element_text(size = base_size, angle = 90),
      axis.text = element_text(size = base_size * 0.9),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 1),
      axis.ticks.length = unit(0.25, "cm"),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.minor.ticks.length = unit(0.15, "cm"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.position = "none",
      plot.title = element_text(size = base_size * 1.1, face = "bold", hjust = 0.5)
    )
}

# ============================================================
# Global Axis Limits
# ============================================================
all_log_abundance <- log10(averaged_intensities %>%
                             select(-UniProt_Entry) %>% unlist() %>% as.numeric())
ecdf_ylim <- range(all_log_abundance, na.rm = TRUE)

# ============================================================
# Generate ECDF Plots and Export as SVG
# ============================================================
for (group in group_definitions) {
  
  group_name <- group$group_name
  mean_col_name <- paste0(group_name, "_mean")
  
  df <- averaged_intensities %>%
    select(UniProt_Entry, abundance = !!sym(mean_col_name)) %>%
    filter(!is.na(abundance) & abundance > 0)
  
  proteins_to_highlight <- highlight_definitions[[group_name]]
  highlights <- df %>%
    filter(UniProt_Entry %in% proteins_to_highlight) %>%
    mutate(protein = ifelse(UniProt_Entry %in% names(protein_names),
                            protein_names[UniProt_Entry],
                            UniProt_Entry))
  
  ecdf_fn <- ecdf(df$abundance)
  highlights <- highlights %>% mutate(percentile = ecdf_fn(abundance))
  ecdf_data <- df %>% arrange(abundance) %>% mutate(percentile = ecdf_fn(abundance))
  
  ecdf_plot <- ggplot(ecdf_data, aes(x = percentile, y = log10(abundance))) +
    geom_line(color = "black", linewidth = 1) +
    geom_point(data = highlights, aes(x = percentile, y = log10(abundance), color = protein),
               size = 5, shape = 16) +
    geom_text_repel(
      data = highlights,
      aes(x = percentile, y = log10(abundance),
          label = paste0(protein, " (", round(percentile*100,1), "%)"),
          color = protein),
      size = 5,
      show.legend = FALSE,
      direction = "both",
      force = 20,
      max.overlaps = Inf,
      segment.color = NA,
      box.padding = 0.5,
      point.padding = 0.5,
      min.segment.length = 0
    ) +
    labs(x = "Percentile", y = expression(log[10]*"(protein abundance intensities)")) +
    scale_color_manual(values = protein_colors) +
    theme_pub() +
    coord_cartesian(ylim = ecdf_ylim)
  
  # Export as SVG
  svg_filename <- paste0("results/figures/ECDF_", group_name, ".svg")
  svglite(svg_filename, width = 7, height = 5)
  print(ecdf_plot)
  dev.off()
}