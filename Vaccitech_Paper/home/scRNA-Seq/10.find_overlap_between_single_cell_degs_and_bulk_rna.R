### THIS FUNCTION:
### 1) FINDS OVERLAP BETWEEN SC DEGs and BULK RNA-SEQ
### GENES THAT ARE UPREGULATED IN BOTH OR DOWNREGULATED IN BOTH ARE "VALIDATED"
### IN THE LARGER COHORT AND ARE PARTICULARLY INTERESTING
### WE ALSO TRACK THE CELL TYPES FROM THE SC DATA
### FOR GENES THAT ARE SIGNIFICANT AT D28 BULK, WE ALSO LOOK AT D2 / D5 / D8 BULK
### WE TRACK WHETHER THESE GENES ARE SIGNIFICANT AT D2 / D5 / D8 AND TRACK FC REGARDLESS
### WE USE ALPHA = 0.05 and ALPHA = 0.1
### 2) WE THEN PLOT SIGNIFICANCE AND FC ACROSS ALL DAYS FOR OUR SET OF GENES THAT PASS
### D28 BULK. WE PLOT D2/D5/D8/D28 AND ALSO D28/D2/D5/D8 (WE ADD LINE BETWEEN 
### D28 AND D2 IN ILLUSTRATOR)
### NOTE: PLOTTING LVL ISN'T USEFUL BECAUSE GENES ARE ONLY SIGNIFICANT AT D28 FOR LVL
### NOTE: PLOTTING MAGICAL GENES ONLY ISN'T USEFUL BECAUSE THERE AREN'T ENOUGH OF THEM
### AND THE FC TREND ISN'T INTERESTING

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Make function that does all of this. Has keyword "innate", "adaptive", and "all"
# This keyword can be used to subset things to innate cell types / adaptive cell types (or just keep all cell types)
# Then just create files with that keyword in directory name
# That way, we can have everything at our fingertips - three simple function calls

# Fill out HVL matrices (upregulated and downregulated genes, with alpha = 0.05 and 0.1)
hvl_upregulated_sc_genes_in_bulk_0.05 <- fill_in_sc_deg_info_for_time_series(innate_sc_pseudobulk_deg_table, high_placebo_counts, high_placebo_metadata,
                                                             paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/"), "up", alpha = 0.05)
saveRDS(hvl_upregulated_sc_genes_in_bulk_0.05, file = paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_0.05.RDS"))

hvl_downregulated_sc_genes_in_bulk_0.05 <- fill_in_sc_deg_info_for_time_series(innate_sc_pseudobulk_deg_table, high_placebo_counts, high_placebo_metadata,
                                                             paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/"), "down", alpha = 0.05)
saveRDS(hvl_downregulated_sc_genes_in_bulk_0.05, file = paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/hvl_downregulated_sc_genes_in_bulk_0.05.RDS"))

#hvl_upregulated_sc_genes_in_bulk_0.1 <- fill_in_sc_deg_info_for_time_series(innate_sc_pseudobulk_deg_table, high_placebo_counts, high_placebo_metadata,
#                                                                     paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/"), "up", alpha = 0.1)
#saveRDS(hvl_upregulated_sc_genes_in_bulk_0.1, file = paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_0.1.RDS"))
#
#hvl_downregulated_sc_genes_in_bulk_0.1 <- fill_in_sc_deg_info_for_time_series(innate_sc_pseudobulk_deg_table, high_placebo_counts, high_placebo_metadata,
#                                                                     paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/"), "down", alpha = 0.1)
#saveRDS(hvl_downregulated_sc_genes_in_bulk_0.1, file = paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/hvl_downregulated_sc_genes_in_bulk_0.1.RDS"))

# Load RDS files
hvl_upregulated_sc_genes_in_bulk_0.05 <- readRDS(paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_0.05.RDS"))
hvl_downregulated_sc_genes_in_bulk_0.05 <- readRDS(paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/hvl_downregulated_sc_genes_in_bulk_0.05.RDS"))
#hvl_upregulated_sc_genes_in_bulk_0.1 <- readRDS(paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_0.1.RDS"))
#hvl_downregulated_sc_genes_in_bulk_0.1 <- readRDS(paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/hvl_downregulated_sc_genes_in_bulk_0.1.RDS"))

# HVL - upregulated, 0.05 alpha
write.table(hvl_upregulated_sc_genes_in_bulk_0.05, file = paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_alpha_0.05.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
high_passing_upregulated_genes_0.05 <- unique(hvl_upregulated_sc_genes_in_bulk_0.05$Gene)

# Plot with D2/D5/D8/D28 in order
hvl_upregulated_sc_genes_in_bulk_0.05_plot <- ggplot(data = hvl_upregulated_sc_genes_in_bulk_0.05, aes(x = Day, y = Gene, size = Fold.Change.Abs, color = Fold.Change.Direction)) +
  geom_point() +
  scale_color_manual(values = c("Positive" = "#00BFC4", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Fold Change of Upregulated Genes (at Day 28) from Innate Immune Cells Across Course of Infection",
    x = "Day (Post Exposure)",
    y = "Gene (Cell Types)",
    size = "Fold Change (Absolute Value)",
    color = "Fold Change Direction"
  ) +
  guides(color  = guide_legend(order = 1),
         size = guide_legend(order = 2)) +
  theme(plot.title = element_text(hjust = 0.6)) + theme(aspect.ratio = 2/1)

ggsave(filename = paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_alpha_0.05_with_cell_types.tiff"), plot = hvl_upregulated_sc_genes_in_bulk_0.05_plot, device='tiff', width = 12, height = 15)

# Alternatively, plot with D28/D2/D5/D8 (can draw a line between D28 and the rest of the days)
hvl_upregulated_sc_genes_in_bulk_0.05_alt <- hvl_upregulated_sc_genes_in_bulk_0.05
hvl_upregulated_sc_genes_in_bulk_0.05_alt$Day <- factor(hvl_upregulated_sc_genes_in_bulk_0.05_alt$Day, levels = c("Day.28","Day.2","Day.5","Day.8"))
hvl_upregulated_sc_genes_in_bulk_0.05_alt_plot <- ggplot(data = hvl_upregulated_sc_genes_in_bulk_0.05_alt, aes(x = Day, y = Gene, size = Fold.Change.Abs, color = Fold.Change.Direction)) +
  geom_point() +
  scale_color_manual(values = c("Positive" = "#00BFC4", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Fold Change of Upregulated Genes (at Day 28) from Innate Immune Cells Across Course of Infection",
    x = "Day (Post Exposure)",
    y = "Gene",
    size = "Fold Change (Absolute Value)",
    color = "Fold Change Direction"
  ) +
  guides(color  = guide_legend(order = 1),
         size = guide_legend(order = 2)) +
  theme(plot.title = element_text(hjust = 0.6)) + theme(aspect.ratio = 3/1)

ggsave(filename = paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_alpha_0.05_alt.png"), plot = hvl_upregulated_sc_genes_in_bulk_0.05_alt_plot, device='png', width = 12, height = 15)

# HVL - downregulated, 0.05 alpha
write.table(hvl_downregulated_sc_genes_in_bulk_0.05, file = paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/hvl_downregulated_sc_genes_in_bulk_alpha_0.05.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
high_passing_downregulated_genes_0.05 <- unique(hvl_downregulated_sc_genes_in_bulk_0.05$Gene)

hvl_downregulated_sc_genes_in_bulk_0.05_plot <- ggplot(data = hvl_downregulated_sc_genes_in_bulk_0.05, aes(x = Day, y = Gene, size = Fold.Change.Abs, color = Fold.Change.Direction)) +
  geom_point() +
  scale_color_manual(values = c("Positive" = "#00BFC4", "Negative" = "#F3756D", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Fold Change of Downregulated Genes (at Day 28) from Innate Immune Cells Across Course of Infection",
    x = "Day (Post Exposure)",
    y = "Gene",
    size = "Fold Change (Absolute Value)",
    color = "Fold Change Direction"
  ) +
  guides(color  = guide_legend(order = 1),
         size = guide_legend(order = 2)) +
  theme(plot.title = element_text(hjust = 0.6)) + theme(aspect.ratio = 2/1)

ggsave(filename = paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/hvl_downregulated_sc_genes_in_bulk_alpha_0.05.tiff"), plot = hvl_downregulated_sc_genes_in_bulk_0.05_plot, device='tiff', width = 12, height = 15)

# Alternatively, plot with D28/D2/D5/D8 (can draw a line between D28 and the rest of the days)
hvl_downregulated_sc_genes_in_bulk_0.05_alt <- hvl_downregulated_sc_genes_in_bulk_0.05
hvl_downregulated_sc_genes_in_bulk_0.05_alt$Day <- factor(hvl_downregulated_sc_genes_in_bulk_0.05_alt$Day, levels = c("Day.28","Day.2","Day.5","Day.8"))
hvl_downregulated_sc_genes_in_bulk_0.05_alt_plot <- ggplot(data = hvl_downregulated_sc_genes_in_bulk_0.05_alt, aes(x = Day, y = Gene, size = Fold.Change.Abs, color = Fold.Change.Direction)) +
  geom_point() +
  scale_color_manual(values = c("Positive" = "#00BFC4", "Negative" = "#F3756D", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Fold Change of Downregulated Genes (at Day 28) from Innate Immune Cells Across Course of Infection",
    x = "Day (Post Exposure)",
    y = "Gene",
    size = "Fold Change (Absolute Value)",
    color = "Fold Change Direction"
  ) +
  guides(color  = guide_legend(order = 1),
         size = guide_legend(order = 2)) +
  theme(plot.title = element_text(hjust = 0.6)) + theme(aspect.ratio = 2/1)

ggsave(filename = paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/hvl_downregulated_sc_genes_in_bulk_alpha_0.05_alt.tiff"), plot = hvl_downregulated_sc_genes_in_bulk_0.05_alt_plot, device='tiff', width = 12, height = 15)

# HVL - upregulated, 0.1 alpha
write.table(hvl_upregulated_sc_genes_in_bulk_0.1, file = paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_alpha_0.1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
high_passing_upregulated_genes_0.1 <- unique(hvl_upregulated_sc_genes_in_bulk_0.1$Gene)

hvl_upregulated_sc_genes_in_bulk_0.1_plot <- ggplot(data = hvl_upregulated_sc_genes_in_bulk_0.1, aes(x = Day, y = Gene, size = Fold.Change.Abs, color = Fold.Change.Direction)) +
  geom_point() +
  scale_color_manual(values = c("Positive" = "#00BFC4", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Fold Change of Upregulated Genes (at Day 28) from Innate Immune Cells Across Course of Infection",
    x = "Day (Post Exposure)",
    y = "Gene",
    size = "Fold Change (Absolute Value)",
    color = "Fold Change Direction"
  ) +
  guides(color  = guide_legend(order = 1),
         size = guide_legend(order = 2)) +
  theme(plot.title = element_text(hjust = 0.6)) + theme(aspect.ratio = 2/1)

ggsave(filename = paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_alpha_0.1.tiff"), plot = hvl_upregulated_sc_genes_in_bulk_0.1_plot, device='tiff', width = 12, height = 15)

# Alternatively, plot with D28/D2/D5/D8 (can draw a line between D28 and the rest of the days)
hvl_upregulated_sc_genes_in_bulk_0.1_alt <- hvl_upregulated_sc_genes_in_bulk_0.1
hvl_upregulated_sc_genes_in_bulk_0.1_alt$Day <- factor(hvl_upregulated_sc_genes_in_bulk_0.1_alt$Day, levels = c("Day.28","Day.2","Day.5","Day.8"))
hvl_upregulated_sc_genes_in_bulk_0.1_alt_plot <- ggplot(data = hvl_upregulated_sc_genes_in_bulk_0.1_alt, aes(x = Day, y = Gene, size = Fold.Change.Abs, color = Fold.Change.Direction)) +
  geom_point() +
  scale_color_manual(values = c("Positive" = "#00BFC4", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Fold Change of Upregulated Genes (at Day 28) from Innate Immune Cells Across Course of Infection",
    x = "Day (Post Exposure)",
    y = "Gene",
    size = "Fold Change (Absolute Value)",
    color = "Fold Change Direction"
  ) +
  guides(color  = guide_legend(order = 1),
         size = guide_legend(order = 2)) +
  theme(plot.title = element_text(hjust = 0.6)) + theme(aspect.ratio = 2/1)

ggsave(filename = paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_alpha_0.1_alt.tiff"), plot = hvl_upregulated_sc_genes_in_bulk_0.1_alt_plot, device='tiff', width = 12, height = 15)

# HVL - downregulated, 0.1 alpha
write.table(hvl_downregulated_sc_genes_in_bulk_0.1, file = paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/hvl_downregulated_sc_genes_in_bulk_alpha_0.1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
high_passing_downregulated_genes_0.1 <- unique(hvl_downregulated_sc_genes_in_bulk_0.1$Gene)

hvl_downregulated_sc_genes_in_bulk_0.1_plot <- ggplot(data = hvl_downregulated_sc_genes_in_bulk_0.1, aes(x = Day, y = Gene, size = Fold.Change.Abs, color = Fold.Change.Direction)) +
  geom_point() +
  scale_color_manual(values = c("Positive" = "#00BFC4", "Negative" = "#F3756D", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Fold Change of Downregulated Genes (at Day 28) from Innate Immune Cells Across Course of Infection",
    x = "Day (Post Exposure)",
    y = "Gene",
    size = "Fold Change (Absolute Value)",
    color = "Fold Change Direction"
  ) +
  guides(color  = guide_legend(order = 1),
         size = guide_legend(order = 2)) +
  theme(plot.title = element_text(hjust = 0.6)) + theme(aspect.ratio = 2/1)

ggsave(filename = paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/hvl_downregulated_sc_genes_in_bulk_alpha_0.1.tiff"), plot = hvl_downregulated_sc_genes_in_bulk_0.1_plot, device='tiff', width = 12, height = 15)

# Alternatively, plot with D28/D2/D5/D8 (can draw a line between D28 and the rest of the days)
hvl_downregulated_sc_genes_in_bulk_0.1_alt <- hvl_downregulated_sc_genes_in_bulk_0.1
hvl_downregulated_sc_genes_in_bulk_0.1_alt$Day <- factor(hvl_downregulated_sc_genes_in_bulk_0.1_alt$Day, levels = c("Day.28","Day.2","Day.5","Day.8"))
hvl_downregulated_sc_genes_in_bulk_0.1_alt_plot <- ggplot(data = hvl_downregulated_sc_genes_in_bulk_0.1_alt, aes(x = Day, y = Gene, size = Fold.Change.Abs, color = Fold.Change.Direction)) +
  geom_point() +
  scale_color_manual(values = c("Positive" = "#00BFC4", "Negative" = "#F3756D", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Fold Change of Downregulated Genes (at Day 28) from Innate Immune Cells Across Course of Infection",
    x = "Day (Post Exposure)",
    y = "Gene",
    size = "Fold Change (Absolute Value)",
    color = "Fold Change Direction"
  ) +
  guides(color  = guide_legend(order = 1),
         size = guide_legend(order = 2)) +
  theme(plot.title = element_text(hjust = 0.6)) + theme(aspect.ratio = 2/1)

ggsave(filename = paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/hvl_downregulated_sc_genes_in_bulk_alpha_0.1_alt.tiff"), plot = hvl_downregulated_sc_genes_in_bulk_0.1_alt_plot, device='tiff', width = 12, height = 15)

# Subset to validated genes
high_all_passing_genes_0.05 <- c(high_passing_upregulated_genes_0.05, high_passing_downregulated_genes_0.05)
innate_sc_pseudobulk_deg_table_passing_hvl_0.05 <- innate_sc_pseudobulk_deg_table[innate_sc_pseudobulk_deg_table$Gene_Name %in% high_all_passing_genes_0.05,]

high_all_passing_genes_0.1 <- c(high_passing_upregulated_genes_0.1, high_passing_downregulated_genes_0.1)
innate_sc_pseudobulk_deg_table_passing_hvl_0.1 <- innate_sc_pseudobulk_deg_table[innate_sc_pseudobulk_deg_table$Gene_Name %in% high_all_passing_genes_0.1,]

# Check validated genes on LVL individuals
lvl_upregulated_sc_genes_in_bulk_0.05 <- fill_in_sc_deg_info_for_time_series(innate_sc_pseudobulk_deg_table_passing_hvl_0.05, low_placebo_counts, low_placebo_metadata,
                                                                             paste0(bulk_results_dir, "lvl_upregulated_sc_genes_found_in_bulk/"), "up", alpha = 0.05)
saveRDS(lvl_upregulated_sc_genes_in_bulk_0.05, file = paste0(bulk_results_dir, "lvl_upregulated_sc_genes_found_in_bulk/lvl_upregulated_sc_genes_in_bulk_0.05.RDS"))

lvl_downregulated_sc_genes_in_bulk_0.05 <- fill_in_sc_deg_info_for_time_series(innate_sc_pseudobulk_deg_table_passing_hvl_0.05, low_placebo_counts, low_placebo_metadata,
                                                                               paste0(bulk_results_dir, "lvl_downregulated_sc_genes_found_in_bulk/"), "down", alpha = 0.05)
saveRDS(lvl_downregulated_sc_genes_in_bulk_0.05, file = paste0(bulk_results_dir, "lvl_downregulated_sc_genes_found_in_bulk/lvl_downregulated_sc_genes_in_bulk_0.05.RDS"))

lvl_upregulated_sc_genes_in_bulk_0.1 <- fill_in_sc_deg_info_for_time_series(innate_sc_pseudobulk_deg_table_passing_hvl_0.1, low_placebo_counts, low_placebo_metadata,
                                                                            paste0(bulk_results_dir, "lvl_upregulated_sc_genes_found_in_bulk/"), "up", alpha = 0.1)
saveRDS(lvl_upregulated_sc_genes_in_bulk_0.1, file = paste0(bulk_results_dir, "lvl_upregulated_sc_genes_found_in_bulk/lvl_upregulated_sc_genes_in_bulk_0.1.RDS"))

lvl_downregulated_sc_genes_in_bulk_0.1 <- fill_in_sc_deg_info_for_time_series(innate_sc_pseudobulk_deg_table_passing_hvl_0.1, low_placebo_counts, low_placebo_metadata,
                                                                              paste0(bulk_results_dir, "lvl_downregulated_sc_genes_found_in_bulk/"), "down", alpha = 0.1)
saveRDS(lvl_downregulated_sc_genes_in_bulk_0.1, file = paste0(bulk_results_dir, "lvl_downregulated_sc_genes_found_in_bulk/lvl_downregulated_sc_genes_in_bulk_0.1.RDS"))

# Load RDS files
lvl_upregulated_sc_genes_in_bulk_0.05 <- readRDS(paste0(bulk_results_dir, "lvl_upregulated_sc_genes_found_in_bulk/lvl_upregulated_sc_genes_in_bulk_0.05.RDS"))
lvl_downregulated_sc_genes_in_bulk_0.05 <- readRDS(paste0(bulk_results_dir, "lvl_downregulated_sc_genes_found_in_bulk/lvl_downregulated_sc_genes_in_bulk_0.05.RDS"))
lvl_upregulated_sc_genes_in_bulk_0.1 <- readRDS(paste0(bulk_results_dir, "lvl_upregulated_sc_genes_found_in_bulk/lvl_upregulated_sc_genes_in_bulk_0.1.RDS"))
lvl_downregulated_sc_genes_in_bulk_0.1 <- readRDS(paste0(bulk_results_dir, "lvl_downregulated_sc_genes_found_in_bulk/lvl_downregulated_sc_genes_in_bulk_0.1.RDS"))

# LVL - upregulated, 0.05 alpha
write.table(lvl_upregulated_sc_genes_in_bulk_0.05, file = paste0(bulk_results_dir, "lvl_upregulated_sc_genes_found_in_bulk/lvl_upregulated_sc_genes_in_bulk_alpha_0.05.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
low_passing_upregulated_genes_0.05 <- unique(lvl_upregulated_sc_genes_in_bulk_0.05$Gene)

# Plot with D2/D5/D8/D28 in order
lvl_upregulated_sc_genes_in_bulk_0.05_plot <- ggplot(data = lvl_upregulated_sc_genes_in_bulk_0.05, aes(x = Day, y = Gene, size = Fold.Change.Abs, color = Fold.Change.Direction)) +
  geom_point() +
  scale_color_manual(values = c("Positive" = "#00BFC4", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Fold Change of Upregulated Genes (at Day 28) from Innate Immune Cells Across Course of Infection for Low Viral Load Individuals",
    x = "Day (Post Exposure)",
    y = "Gene",
    size = "Fold Change (Absolute Value)",
    color = "Fold Change Direction"
  ) +
  guides(color  = guide_legend(order = 1),
         size = guide_legend(order = 2)) +
  theme(plot.title = element_text(hjust = 0.6)) + theme(aspect.ratio = 2/1)

ggsave(filename = paste0(bulk_results_dir, "lvl_upregulated_sc_genes_found_in_bulk/lvl_upregulated_sc_genes_in_bulk_alpha_0.05.tiff"), plot = lvl_upregulated_sc_genes_in_bulk_0.05_plot, device='tiff', width = 12, height = 15)

# Current thought: Use plots for HVL individuals (up and down regulated). 
# Don't use plots for LVL individuals because only D28 is significant. Instead, just add the LVL D28 significance as a side note.

# Possible idea: Try subset of genes that are found in MAGICAL circuits and plot those
# Upregulated, alpha = 0.05. Doesn't really work - too few genes.
magical_output_dir <- paste0(sc_magical_dir, "Output/")
overall_magical_results <- read.table(paste0(magical_output_dir, "MAGICAL_overall_output.tsv"), sep = "\t", header = TRUE)
innate_cell_types_with_underscore <- gsub(" ", "_", innate_cell_types)
innate_magical_results <- overall_magical_results[overall_magical_results$Cell_Type %in% innate_cell_types_with_underscore,]
innate_magical_genes <- unique(innate_magical_results$Gene_symbol)
hvl_upregulated_sc_genes_in_bulk_0.05_magical <- hvl_upregulated_sc_genes_in_bulk_0.05[hvl_upregulated_sc_genes_in_bulk_0.05$Gene %in% innate_magical_genes,]

# Upregulated, alpha = 0.1. Eh, maybe a little better but not much - too few genes.
hvl_upregulated_sc_genes_in_bulk_0.1_magical <- hvl_upregulated_sc_genes_in_bulk_0.1[hvl_upregulated_sc_genes_in_bulk_0.1$Gene %in% innate_magical_genes,]

# Fill out HVL matrices for ALL genes (combined cell types) for motif overlap (script 4)
hvl_upregulated_sc_genes_in_bulk_combined_cell_types_0.05 <- fill_in_sc_deg_info_for_time_series(sc_pseudobulk_deg_combined_cell_types_table, high_placebo_counts, high_placebo_metadata,
                                                                             paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk_combined_cell_types/"), "up", alpha = 0.05)
saveRDS(hvl_upregulated_sc_genes_in_bulk_combined_cell_types_0.05, file = paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk_combined_cell_types/hvl_upregulated_sc_genes_in_bulk_combined_cell_types_0.05.RDS"))

hvl_downregulated_sc_genes_in_bulk_combined_cell_types_0.05 <- fill_in_sc_deg_info_for_time_series(sc_pseudobulk_deg_combined_cell_types_table, high_placebo_counts, high_placebo_metadata,
                                                                               paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk_combined_cell_types/"), "down", alpha = 0.05)
saveRDS(hvl_downregulated_sc_genes_in_bulk_combined_cell_types_0.05, file = paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk_combined_cell_types/hvl_downregulated_sc_genes_in_bulk_combined_cell_types_0.05.RDS"))

# Load HVL matrices for ALL genes (combined cell types)
hvl_upregulated_sc_genes_in_bulk_combined_cell_types_0.05 <- readRDS(paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk_combined_cell_types/hvl_upregulated_sc_genes_in_bulk_combined_cell_types_0.05.RDS"))
hvl_downregulated_sc_genes_in_bulk_combined_cell_types_0.05 <- readRDS(paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk_combined_cell_types/hvl_downregulated_sc_genes_in_bulk_combined_cell_types_0.05.RDS"))



#### here is some code to plot adjusted p-value with color gradient - not sure I like it better, and don't know
#### how to handle negative and positive fold change

library(dplyr)

# Convert Adjusted.P.Value to numeric (if it's not already numeric)
hvl_upregulated_sc_genes_in_bulk_0.05 <- hvl_upregulated_sc_genes_in_bulk_0.05 %>%
  mutate(Adjusted.P.Value = as.numeric(Adjusted.P.Value))

# Plot with color gradient based on Adjusted.P.Value for positive Fold Change
ggplot(data = hvl_upregulated_sc_genes_in_bulk_0.05, aes(x = Day, y = Gene, size = Fold.Change.Abs)) +
  geom_point(aes(color = ifelse(Fold.Change.Direction == "Positive", Adjusted.P.Value, NA)), show.legend = TRUE) +
  scale_color_gradientn(colors = c("aquamarine4", "aquamarine"), na.value = "grey") +
  theme_minimal() +
  labs(
    title = "Fold Change of Upregulated Genes (at Day 28) from Innate Immune Cells Across Course of Infection",
    x = "Day (Post Exposure)",
    y = "Gene (Cell Types)",
    size = "Fold Change (Absolute Value)",
    color = "Adjusted P-Value (Significant)"
  ) +
  guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2)) +
  theme(plot.title = element_text(hjust = 0.6)) +
  theme(aspect.ratio = 2/1)

