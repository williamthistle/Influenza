# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# This method finds overlap between enriched transcription factors in mintchip and those genes expressed in the scRNA-seq data
# We also flag whether the gene was also found to be significant in bulk
overlapping_motifs_mintchip <- find_overlapping_motifs_between_mintchip_and_rna(homer_dir, sc_pseudobulk_deg_combined_cell_types_table, possible_markers, unique(hvl_upregulated_sc_genes_in_bulk_combined_cell_types_0.05$Gene), unique(hvl_downregulated_sc_genes_in_bulk_combined_cell_types_0.05$Gene))
# Only one TF (ETS1) found in bulk, but a few were found in single cell data. Also interesting!