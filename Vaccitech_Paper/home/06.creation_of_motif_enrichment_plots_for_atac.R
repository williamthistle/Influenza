# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

atac_tf_motif_enrichment_dir <- paste0(sc_das_dir, "diff_peaks/")
atac_innate_cell_types <- c("CD14_Mono", "CD16_Mono", "NK")

atac_tf_motif_enrichment_plot_dir <- paste0(sc_das_dir, "motif_enrichment_plots/")
if (!dir.exists(atac_tf_motif_enrichment_plot_dir)) {dir.create(atac_tf_motif_enrichment_plot_dir)}

# Cell Type
# Motif
# -log 10 p-value
# First, do upregulated motifs - top 15
up_atac_tf_motif_enrichment_df <- data.frame(TF = character(), Enrichment = numeric(), Negative.Log.Adjusted.P.Value <- numeric(), Cell.Type = character())
for(atac_innate_cell_type in atac_innate_cell_types) {
  current_tf_motif_enrichment_results <- read.table(paste0(atac_tf_motif_enrichment_dir, atac_innate_cell_type, "_D28_D1_motif_up_pseudobulk_only.tsv"), sep = "\t", header = TRUE)
  current_tf_motif_enrichment_results$TF <- sub("_.*", "", current_tf_motif_enrichment_results$TF)
  current_tf_motif_enrichment_results$Cell_Type <- sub("_", " ", atac_innate_cell_type)
  current_tf_motif_enrichment_results_subset <- current_tf_motif_enrichment_results[1:15,]
  current_tf_motif_enrichment_results_subset <- current_tf_motif_enrichment_results_subset[,-c(4, 5)]
  colnames(current_tf_motif_enrichment_results_subset) <- c("TF", "Enrichment", "Negative.Log.Adjusted.P.Value", "Cell.Type")
  up_atac_tf_motif_enrichment_df <- rbind(up_atac_tf_motif_enrichment_df, current_tf_motif_enrichment_results_subset)
}

up_atac_tf_motif_enrichment_df_plot <- ggplot(data = up_atac_tf_motif_enrichment_df, aes(x = Cell.Type, y = TF, size = Enrichment, color = Negative.Log.Adjusted.P.Value)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Fold Change of Upregulated Genes (at Day 28) from Innate Immune Cells Across Course of Infection",
    x = "Day (Post Exposure)",
    y = "Gene",
    size = "Enrichment",
    color = "-Log(Adjusted P Value)"
  ) +
  theme(plot.title = element_text(hjust = 0.6)) + theme(aspect.ratio = 2/1)

ggsave(filename = paste0(atac_tf_motif_enrichment_plot_dir, "up_atac_tf_motif_enrichment_df_plot.tiff"), plot = up_atac_tf_motif_enrichment_df_plot, device='tiff', dpi = 300)

# Second, do downregulated motifs - top 15
down_atac_tf_motif_enrichment_df <- data.frame(TF = character(), Enrichment = numeric(), Negative.Log.Adjusted.P.Value <- numeric(), Cell.Type = character())
for(atac_innate_cell_type in atac_innate_cell_types) {
  current_tf_motif_enrichment_results <- read.table(paste0(atac_tf_motif_enrichment_dir, atac_innate_cell_type, "_D28_D1_motif_down_pseudobulk_only.tsv"), sep = "\t", header = TRUE)
  current_tf_motif_enrichment_results$TF <- sub("_.*", "", current_tf_motif_enrichment_results$TF)
  current_tf_motif_enrichment_results$Cell_Type <- sub("_", " ", atac_innate_cell_type)
  current_tf_motif_enrichment_results_subset <- current_tf_motif_enrichment_results[1:15,]
  current_tf_motif_enrichment_results_subset <- current_tf_motif_enrichment_results_subset[,-c(4, 5)]
  colnames(current_tf_motif_enrichment_results_subset) <- c("TF", "Enrichment", "Negative.Log.Adjusted.P.Value", "Cell.Type")
  down_atac_tf_motif_enrichment_df <- rbind(down_atac_tf_motif_enrichment_df, current_tf_motif_enrichment_results_subset)
}

down_atac_tf_motif_enrichment_df_plot <- ggplot(data = down_atac_tf_motif_enrichment_df, aes(x = Cell.Type, y = TF, size = Enrichment, color = Negative.Log.Adjusted.P.Value)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Fold Change of downregulated Genes (at Day 28) from Innate Immune Cells Across Course of Infection",
    x = "Day (Post Exposure)",
    y = "Gene",
    size = "Enrichment",
    color = "-Log(Adjusted P Value)"
  ) +
  theme(plot.title = element_text(hjust = 0.6)) + theme(aspect.ratio = 2/1)

ggsave(filename = paste0(atac_tf_motif_enrichment_plot_dir, "down_atac_tf_motif_enrichment_df_plot.tiff"), plot = down_atac_tf_motif_enrichment_df_plot, device='tiff', dpi = 300)
