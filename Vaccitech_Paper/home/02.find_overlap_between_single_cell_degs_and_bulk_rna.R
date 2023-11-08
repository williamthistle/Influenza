# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Fill out HVL matrices (upregulated and downregulated genes, with alpha = 0.05 and 0.1)
hvl_upregulated_sc_genes_in_bulk_0.05 <- fill_in_sc_deg_info_for_time_series(innate_sc_pseudobulk_deg_table, high_placebo_counts, high_placebo_metadata,
                                                             paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/"), "up", alpha = 0.05)
saveRDS(hvl_upregulated_sc_genes_in_bulk_0.05, file = paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_0.05.RDS"))
hvl_downregulated_sc_genes_in_bulk_0.05 <- fill_in_sc_deg_info_for_time_series(innate_sc_pseudobulk_deg_table, high_placebo_counts, high_placebo_metadata,
                                                             paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/"), "down", alpha = 0.05)
saveRDS(hvl_downregulated_sc_genes_in_bulk_0.05, file = paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/hvl_downregulated_sc_genes_in_bulk_0.05.RDS"))
hvl_upregulated_sc_genes_in_bulk_0.1 <- fill_in_sc_deg_info_for_time_series(innate_sc_pseudobulk_deg_table, high_placebo_counts, high_placebo_metadata,
                                                                     paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/"), "up", alpha = 0.1)
saveRDS(hvl_upregulated_sc_genes_in_bulk_0.1, file = paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_0.1.RDS"))
hvl_downregulated_sc_genes_in_bulk_0.1 <- fill_in_sc_deg_info_for_time_series(innate_sc_pseudobulk_deg_table, high_placebo_counts, high_placebo_metadata,
                                                                     paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/"), "down", alpha = 0.1)
saveRDS(hvl_downregulated_sc_genes_in_bulk_0.1, file = paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/hvl_downregulated_sc_genes_in_bulk_0.1.RDS"))

# HVL - upregulated, 0.05 alpha
write.table(hvl_upregulated_sc_genes_in_bulk_0.05, file = paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_alpha_0.05.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
high_passing_upregulated_genes <- unique(hvl_upregulated_sc_genes_in_bulk_0.05$Gene)

# Plot with D2/D5/D8/D28 in order
hvl_upregulated_sc_genes_in_bulk_0.05_plot <- ggplot(data = hvl_upregulated_sc_genes_in_bulk_0.05, aes(x = Day, y = Gene, size = Fold.Change.Abs, color = Fold.Change.Direction)) +
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

ggsave(filename = paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_alpha_0.05.tiff"), plot = hvl_upregulated_sc_genes_in_bulk_0.05_plot, device='tiff', dpi=300)

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
  theme(plot.title = element_text(hjust = 0.6)) + theme(aspect.ratio = 2/1)

ggsave(filename = paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_alpha_0.05_alt.tiff"), plot = hvl_upregulated_sc_genes_in_bulk_0.05_alt_plot, device='tiff', width = 12, height = 15)

# HVL - downregulated, 0.05 alpha
write.table(hvl_downregulated_sc_genes_in_bulk_0.05, file = paste0(bulk_results_dir, "hvl_downregulated_sc_genes_found_in_bulk/hvl_downregulated_sc_genes_in_bulk_alpha_0.05.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
high_passing_downregulated_genes <- unique(hvl_downregulated_sc_genes_in_bulk_0.05$Gene)

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
high_passing_upregulated_genes <- unique(hvl_upregulated_sc_genes_in_bulk_0.1$Gene)

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
high_passing_downregulated_genes <- unique(hvl_downregulated_sc_genes_in_bulk_0.1$Gene)

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












passing_bulk_genes <- c(high_passing_pos_genes, high_passing_neg_genes)

# Check these validated genes on low viral load individuals
low_pos_pseudobulk_sc_DEGs_bulk_passing_df <- find_degs_across_time_points_for_gene_list(raw_low_placebo_period_2_D2_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D5_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D8_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D28_vs_D_minus_1_results,
                                                                                          low_pos_pseudobulk_sc_DEGs_bulk_passing_df, sc_pseudobulk_DEG_table[sc_pseudobulk_DEG_table$Gene_Name %in% high_passing_pos_genes,])
low_neg_pseudobulk_sc_DEGs_bulk_passing_df <- find_degs_across_time_points_for_gene_list(raw_low_placebo_period_2_D2_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D5_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D8_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D28_vs_D_minus_1_results,
                                                                                          low_neg_pseudobulk_sc_DEGs_bulk_passing_df, sc_pseudobulk_DEG_table[sc_pseudobulk_DEG_table$Gene_Name %in% high_passing_neg_genes,])
# pos: 6 genes
low_pos_pseudobulk_sc_DEGs_bulk_passing_df <- fill_in_special_notes(low_pos_pseudobulk_sc_DEGs_bulk_passing_df)
write.table(low_pos_pseudobulk_sc_DEGs_bulk_passing_df, file = paste0(onedrive_dir, "Vaccitech_Paper/low_passing_pos_df.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
low_passing_pos_genes <- low_pos_pseudobulk_sc_DEGs_bulk_passing_df[low_pos_pseudobulk_sc_DEGs_bulk_passing_df$D28_0.2 == TRUE,]$gene
# neg: 19 genes
low_neg_pseudobulk_sc_DEGs_bulk_passing_df <- fill_in_special_notes(low_neg_pseudobulk_sc_DEGs_bulk_passing_df)
write.table(low_neg_pseudobulk_sc_DEGs_bulk_passing_df, file = paste0(onedrive_dir, "Vaccitech_Paper/low_passing_neg_df.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
low_passing_neg_genes <- low_neg_pseudobulk_sc_DEGs_bulk_passing_df[low_neg_pseudobulk_sc_DEGs_bulk_passing_df$D28_negative_0.2 == TRUE,]$gene



# FLAGGED GENES - those genes which have negative FC in SC data but positive FC in bulk data (or vice versa)
# Not relevant with new formatting, but can maybe replicate it with small changes
#flagged_high_pos_genes <- high_pos_pseudobulk_sc_DEGs_bulk_passing_df[high_pos_pseudobulk_sc_DEGs_bulk_passing_df$D28_negative_0.2 == TRUE,]$gene
#flagged_high_pos_gene_df <- sc_pseudobulk_DEG_table[sc_pseudobulk_DEG_table$Gene_Name %in% flagged_high_pos_genes,]

#flagged_high_neg_genes <- high_neg_pseudobulk_sc_DEGs_bulk_passing_df[high_neg_pseudobulk_sc_DEGs_bulk_passing_df$D28_0.2 == TRUE,]$gene
#flagged_high_neg_gene_df <- sc_pseudobulk_DEG_table[sc_pseudobulk_DEG_table$Gene_Name %in% flagged_high_neg_genes,]

#flagged_low_pos_genes <- low_pos_pseudobulk_sc_DEGs_bulk_passing_df[low_pos_pseudobulk_sc_DEGs_bulk_passing_df$D28_negative_0.2 == TRUE,]$gene
#flagged_low_pos_gene_df <- sc_pseudobulk_DEG_table[sc_pseudobulk_DEG_table$Gene_Name %in% flagged_low_pos_genes,]

#flagged_low_neg_genes <- low_neg_pseudobulk_sc_DEGs_bulk_passing_df[low_neg_pseudobulk_sc_DEGs_bulk_passing_df$D28_0.2 == TRUE,]$gene
#flagged_low_neg_gene_df <- sc_pseudobulk_DEG_table[sc_pseudobulk_DEG_table$Gene_Name %in% flagged_low_neg_genes,]




