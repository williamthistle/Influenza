# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

possible_cell_types <- c("CD4_Naive", "CD8_Naive", "CD4_Memory", "CD8_Memory", "cDC", "HSPC", "pDC", "Platelet", "Plasmablast", "Proliferating", "NK", "NK_CD56bright", "T_Naive", "CD14_Mono", "CD16_Mono", "MAIT", "B", "B_naive", "B_memory")

# This method finds overlap between enriched transcription factors in pseudobulk ATAC-seq and those genes expressed in the scRNA-seq data
# We also flag whether the gene was also found to be significant in bulk
overlapping_motifs_sc <- find_overlapping_motifs_between_atac_and_rna(paste0(sc_das_dir, "diff_peaks/"), sc_pseudobulk_deg_combined_cell_types_table, possible_cell_types, unique(hvl_upregulated_sc_genes_in_bulk_combined_cell_types_0.05$Gene), unique(hvl_downregulated_sc_genes_in_bulk_combined_cell_types_0.05$Gene))
# None of the enriched TFs were found in bulk, but a few were found in the single cell data. This is interesting, I think!