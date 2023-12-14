# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

atac_tf_motif_enrichment_dir <- paste0(sc_das_dir, "motif_enrichment/")
atac_innate_cell_types <- c("CD14_Mono", "CD16_Mono")

atac_tf_motif_enrichment_plot_dir <- paste0(sc_das_dir, "motif_enrichment_plots/")
if (!dir.exists(atac_tf_motif_enrichment_plot_dir)) {dir.create(atac_tf_motif_enrichment_plot_dir)}

# Cell Type
# Motif
# -log 10 p-value
# First, do upregulated motifs - top 15
filtered_rows <- list()
for(atac_innate_cell_type in atac_innate_cell_types) {
  if(atac_innate_cell_type == "CD14_Mono") {
    peak_count <- "total_peaks_2166"
  } else {
    peak_count <- "total_peaks_3704"
  }
  current_tf_motif_enrichment_results <- read.table(paste0(atac_tf_motif_enrichment_dir, atac_innate_cell_type, "/sc_filtered/0.01/with_bg/D28-vs-D_minus_1-degs-", atac_innate_cell_type, "-sc_filtered_pct_0.01_FC_1_", peak_count, "_pos_motifs_with_bg.tsv"), sep = "\t", header = TRUE)
  current_tf_motif_enrichment_results$Cell.Type <- sub("_", " ", atac_innate_cell_type)
  current_tf_motif_enrichment_results_subset <- current_tf_motif_enrichment_results[1:15,]
  # current_tf_motif_enrichment_results_subset <- current_tf_motif_enrichment_results_subset[,-c(4, 5)]
  # colnames(current_tf_motif_enrichment_results_subset) <- c("TF", "Enrichment", "Negative.Log.Adjusted.P.Value", "Cell.Type")
  # up_atac_tf_motif_enrichment_df <- rbind(up_atac_tf_motif_enrichment_df, current_tf_motif_enrichment_results_subset)
  filtered_rows[[length(filtered_rows) + 1]] <- current_tf_motif_enrichment_results_subset
}

up_atac_tf_motif_enrichment_df <- do.call(rbind, filtered_rows)
up_atac_tf_motif_enrichment_df$Negative.Log.Adjusted.P.Value <- -log10(up_atac_tf_motif_enrichment_df$p.adjust)

up_atac_tf_motif_enrichment_df_plot <- ggplot(data = up_atac_tf_motif_enrichment_df, aes(x = Cell.Type, y = motif.name, size = fold.enrichment, color = Negative.Log.Adjusted.P.Value)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "TF Motif Enrichment for Post-Exposure Sites in Monocytes",
    x = "Cell Type",
    y = "TF",
    size = "Enrichment",
    color = "-Log(Adjusted P Value)"
  ) +
  theme(plot.title = element_text(hjust = 0.6)) + theme(aspect.ratio = 2/1)

ggsave(filename = paste0(atac_tf_motif_enrichment_plot_dir, "up_atac_tf_motif_enrichment_df_plot.tiff"), plot = up_atac_tf_motif_enrichment_df_plot, device='tiff', dpi = 300)

# Second, do downregulated motifs - top 15
filtered_rows_negative <- list()
for(atac_innate_cell_type in atac_innate_cell_types) {
  if(atac_innate_cell_type == "CD14_Mono") {
    peak_count <- "total_peaks_2397"
  } else {
    peak_count <- "total_peaks_2838"
  }
  current_tf_motif_enrichment_results <- read.table(paste0(atac_tf_motif_enrichment_dir, atac_innate_cell_type, "/sc_filtered/0.01/with_bg/D28-vs-D_minus_1-degs-", atac_innate_cell_type, "-sc_filtered_pct_0.01_FC_-1_", peak_count, "_neg_motifs_with_bg.tsv"), sep = "\t", header = TRUE)
  current_tf_motif_enrichment_results$Cell.Type <- sub("_", " ", atac_innate_cell_type)
  current_tf_motif_enrichment_results_subset <- current_tf_motif_enrichment_results[1:15,]
  # current_tf_motif_enrichment_results_subset <- current_tf_motif_enrichment_results_subset[,-c(4, 5)]
  # colnames(current_tf_motif_enrichment_results_subset) <- c("TF", "Enrichment", "Negative.Log.Adjusted.P.Value", "Cell.Type")
  # up_atac_tf_motif_enrichment_df <- rbind(up_atac_tf_motif_enrichment_df, current_tf_motif_enrichment_results_subset)
  filtered_rows_negative[[length(filtered_rows_negative) + 1]] <- current_tf_motif_enrichment_results_subset
}

down_atac_tf_motif_enrichment_df <- do.call(rbind, filtered_rows_negative)
down_atac_tf_motif_enrichment_df$Negative.Log.Adjusted.P.Value <- -log10(down_atac_tf_motif_enrichment_df$p.adjust)

down_atac_tf_motif_enrichment_df_plot <- ggplot(data = down_atac_tf_motif_enrichment_df, aes(x = Cell.Type, y = motif.name, size = fold.enrichment, color = Negative.Log.Adjusted.P.Value)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "TF Motif Enrichment for Differentially Accessible Sites in Monocytes (Pre-Exposure)",
    x = "Cell Type",
    y = "TF",
    size = "Enrichment",
    color = "-Log(Adjusted P Value)"
  ) +
  theme(plot.title = element_text(hjust = 0.6)) + theme(aspect.ratio = 2/1)


ggsave(filename = paste0(atac_tf_motif_enrichment_plot_dir, "down_atac_tf_motif_enrichment_df_plot.tiff"), plot = down_atac_tf_motif_enrichment_df_plot, device='tiff', dpi = 300)
