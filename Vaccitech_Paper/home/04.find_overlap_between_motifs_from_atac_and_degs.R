# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

overlapping_motifs_sc <- find_overlapping_motifs_between_atac_and_rna(paste0(sc_peak_dir, "diff_peaks/"), sc_pseudobulk_deg_combined_cell_types_table, possible_cell_types, high_pos_pseudobulk_sc_genes_bulk_passing_df$gene, high_neg_pseudobulk_sc_genes_bulk_passing_df$gene)
