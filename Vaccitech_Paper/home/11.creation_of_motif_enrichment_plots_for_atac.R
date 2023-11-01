# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

atac_tf_motif_enrichment_dir <- paste0(sc_das_dir, "diff_peaks/")
atac_innate_cell_types <- c("CD14_Mono", "CD16_Mono", "NK")

# Cell Type
# Motif
# -log 10 p-value
# First, do upregulated motifs - top 10?
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

