# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# STEP 1: SET UP MAGICAL INPUT FILES BASED ON RNA / ATAC DATA

# a) Cell type candidate genes.txt
cell_type_candidate_gene_dir <- paste0(sc_magical_dir, "Candidate_Genes/")
for(cell_type in unique(sc_pseudobulk_deg_combined_cell_types_table$Cell_Type)) {
  write.table(sc_pseudobulk_deg_combined_cell_types_table[sc_pseudobulk_deg_combined_cell_types_table$Cell_Type == cell_type,]$Gene_Name,
              file = paste0(cell_type_candidate_gene_dir, sub(" ", "_", cell_type), "_Candidate_Genes.txt"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

# b) Cell type candidate peaks.txt
cell_type_candidate_peak_dir <- paste0(sc_magical_dir, "Candidate_Peaks/")
for(cell_type in unique(sc_das_lenient$Cell_Type)) {
  cell_type_sc_das_lenient <- sc_das_lenient[sc_das_lenient$Cell_Type == cell_type,]
  cell_type_sc_das_lenient <- cell_type_sc_das_lenient[,c("chr", "start", "end")]
  write.table(cell_type_sc_das_lenient,
              file = paste0(cell_type_candidate_peak_dir, sub(" ", "_", cell_type), "_Candidate_Peaks.txt"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

# c) Cell type scATAC cell meta.txt
cell_type_scATAC_cell_metadata_dir <- paste0(sc_magical_dir, "scATAC_Cell_Metadata/")
for(cell_type in unique(atac_cell_metadata$cell_type)) {
  cell_type_atac_cell_metadata <- atac_cell_metadata[atac_cell_metadata$cell_type == cell_type,]
  write.table(cell_type_atac_cell_metadata,
              file = paste0(cell_type_scATAC_cell_metadata_dir, sub(" ", "_", cell_type), "_scATAC_Cell_Metadata.txt"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

# d) Cell type scATAC read count.txt


# e) Cell type scRNA cell meta.txt
# f) Cell type scRNA read count.txt
# g) hg38_Refseq.txt
hg38_ref_seq <- paste0(sc_magical_dir, "hg38_Refseq.txt")
# h) Motif mapping prior.txt
# i) Motifs.txt
motifs <- paste0(sc_magical_dir, "Motifs.txt")
# j) RaoGM12878_40kb_TopDomTADs_filtered_hg38.txt
tads <- paste0(sc_magical_dir, "RaoGM12878_40kb_TopDomTADs_filtered_hg38.txt")
# k) scATAC peaks.txt
# l) scRNA genes.txt
