#hvl_upregulated_sc_genes_in_bulk_unfiltered <- fill_in_sc_deg_info_for_time_series(innate_sc_pseudobulk_deg_table, high_placebo_counts, high_placebo_metadata,
#                                                                             paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk_unfiltered/"), "up", alpha = 0.05,
#                                                                             filter_D28 = FALSE)
#saveRDS(hvl_upregulated_sc_genes_in_bulk_unfiltered, file = paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_unfiltered.RDS"))

hvl_upregulated_sc_genes_in_bulk_unfiltered <- readRDS(paste0(bulk_results_dir, "hvl_upregulated_sc_genes_found_in_bulk/hvl_upregulated_sc_genes_in_bulk_unfiltered.RDS"))

# Remove cell types
hvl_upregulated_sc_genes_in_bulk_unfiltered <- hvl_upregulated_sc_genes_in_bulk_unfiltered %>%
  mutate(Gene = sapply(Gene, remove_text_in_parentheses))

# Change Positive / Negative to Significant
hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction <- replace(hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction,
                                                                       hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction == "Positive",
                                                                       "Significant")

hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction <- replace(hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction,
                                                                       hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction == "Negative",
                                                                       "Significant")

# Change absolute fold change to include direction of fold change (and rename column to remove Abs)
hvl_upregulated_sc_genes_in_bulk_unfiltered <- hvl_upregulated_sc_genes_in_bulk_unfiltered %>%
  mutate(Fold.Change = ifelse(Fold.Change.Direction.Raw == "Negative", -Fold.Change.Abs, Fold.Change.Abs)) %>%
  select(-Fold.Change.Abs)

hvl_upregulated_sc_genes_in_bulk_unfiltered$Day <- factor(hvl_upregulated_sc_genes_in_bulk_unfiltered$Day, levels = c("Day.2", "Day.5", "Day.8", "Day.28"))
hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction.Raw <- factor(hvl_upregulated_sc_genes_in_bulk_unfiltered$Fold.Change.Direction.Raw, 
                                                                                levels = c("Positive", "Negative"))

hvl_upregulated_sc_genes_in_bulk_unfiltered_plot <- ggplot(hvl_upregulated_sc_genes_in_bulk_unfiltered, aes(Day, Fold.Change, color = Fold.Change.Direction.Raw)) +
  geom_jitter() + theme_minimal(base_size = 14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("Day (Post-Exposure)") +
  ylab("log2FC (Gene Expression)") + ggtitle("Fold Change Over Time for Single Cell DEGs Upregulated at Day 28") +
  theme(plot.title = element_text(hjust = 0.5)) + guides(color=guide_legend(title="Fold Change Direction"))

ggsave(filename = paste0("C:/Users/willi/Desktop/", "bulk_deg_jitterplot.tiff"), plot = bulk_degs_jitter_plot, device='tiff', dpi=300)


