# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

source(paste0(base_dir, "extra_functions/MAGICAL_functions.R"))

# STEP 0: Set up directories and global files
cell_type_candidate_gene_dir <- paste0(sc_magical_dir, "Candidate_Genes/")
if (!dir.exists(cell_type_candidate_gene_dir)) {dir.create(cell_type_candidate_gene_dir, recursive = TRUE)}
cell_type_candidate_peak_overlap_dir <- paste0(sc_magical_dir, "Candidate_Peaks/Overlap_SC_Pseudo/")
if (!dir.exists(cell_type_candidate_peak_overlap_dir)) {dir.create(cell_type_candidate_peak_overlap_dir, recursive = TRUE)}
cell_type_candidate_peak_pseudo_dir <- paste0(sc_magical_dir, "Candidate_Peaks/Pseudo/")
if (!dir.exists(cell_type_candidate_peak_pseudo_dir)) {dir.create(cell_type_candidate_peak_pseudo_dir, recursive = TRUE)}
cell_type_scATAC_cell_metadata_dir <- paste0(sc_magical_dir, "scATAC_Cell_Metadata/")
if (!dir.exists(cell_type_scATAC_cell_metadata_dir)) {dir.create(cell_type_scATAC_cell_metadata_dir, recursive = TRUE)}
cell_type_scRNA_cell_metadata_dir <- paste0(sc_magical_dir, "scRNA_Cell_Metadata/")
if (!dir.exists(cell_type_scRNA_cell_metadata_dir)) {dir.create(cell_type_scRNA_cell_metadata_dir, recursive = TRUE)}
cell_type_scATAC_read_counts_dir <- paste0(sc_magical_dir, "scATAC_Read_Counts/")
if (!dir.exists(cell_type_scATAC_read_counts_dir)) {dir.create(cell_type_scATAC_read_counts_dir, recursive = TRUE)}
cell_type_scRNA_read_counts_dir <- paste0(sc_magical_dir, "scRNA_Read_Counts/")
if (!dir.exists(cell_type_scRNA_read_counts_dir)) {dir.create(cell_type_scRNA_read_counts_dir, recursive = TRUE)}
magical_output_dir <- paste0(sc_magical_dir, "Output/")
if (!dir.exists(magical_output_dir)) {dir.create(magical_output_dir, recursive = TRUE)}

hg38_ref_seq <- paste0(sc_magical_dir, "hg38_Refseq.txt")
peak_coordinates <- paste0(sc_magical_dir, "scATAC_Peak_Coordinates/HVL_ATAC_peak_coordinates.tsv")
gene_list <- paste0(sc_magical_dir, "scRNA_Genes/HVL_RNA_genes.tsv")
motifs <- paste0(sc_magical_dir, "Motifs.txt")
motif_prior <- paste0(sc_magical_dir, "scATAC_Motif_Mapping_Prior/HVL_ATAC_motif_mapping_prior.tsv")
tads <- paste0(sc_magical_dir, "RaoGM12878_40kb_TopDomTADs_filtered_hg38.txt")
distance_control=5e5

# STEP 1: CREATE CANDIDATE GENE AND PEAK FILES

# a) Cell type candidate genes.txt

for(cell_type in unique(sc_pseudobulk_deg_combined_cell_types_table$Cell_Type)) {
  write.table(sc_pseudobulk_deg_combined_cell_types_table[sc_pseudobulk_deg_combined_cell_types_table$Cell_Type == cell_type,]$Gene_Name,
              file = paste0(cell_type_candidate_gene_dir, sub(" ", "_", cell_type), "_Candidate_Genes.txt"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

# b) Cell type candidate peaks.txt (overlap)
for(cell_type in unique(sc_das_lenient$Cell_Type)) {
  cell_type_sc_das_lenient <- sc_das_lenient[sc_das_lenient$Cell_Type == cell_type,]
  cell_type_sc_das_lenient <- cell_type_sc_das_lenient[,c("chr", "start", "end")]
  write.table(cell_type_sc_das_lenient,
              file = paste0(cell_type_candidate_peak_overlap_dir, sub(" ", "_", cell_type), "_Candidate_Peaks.txt"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

# c) Cell type candidate peaks.txt (pseudo)
for(cell_type in unique(sc_das_lenient$Cell_Type)) {
  current_pseudobulk_das <- read.table(paste0(sc_das_dir, "diff_peaks/", sub(" ", "_", cell_type), "_D28_D1_diff_pseudo_filtered.tsv"), sep = "\t",
                                       header = TRUE)
  current_pseudobulk_das <- current_pseudobulk_das[,c("chr", "start", "end")]
  write.table(current_pseudobulk_das,
              file = paste0(cell_type_candidate_peak_pseudo_dir, sub(" ", "_", cell_type), "_Candidate_Peaks.txt"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

# STEP 2: FIND MATCHING CELL TYPES FROM RNA / ATAC FOR MAGICAL ANALYSIS
possible_RNA_types <- list.files(path = cell_type_candidate_gene_dir)
possible_RNA_types <- unlist(strsplit(possible_RNA_types, "_Candidate_Genes.txt"))

possible_ATAC_types <- list.files(path = cell_type_candidate_peak_pseudo_dir)
possible_ATAC_types <- unlist(strsplit(possible_ATAC_types, "_Candidate_Peaks.txt"))

cell_types <- intersect(possible_RNA_types, possible_ATAC_types)

# STEP 3: RUN MAGICAL
set.seed(get_speedi_seed())
for(cell_type in cell_types) {
  print(cell_type)
  current_candidate_genes <- paste0(cell_type_candidate_gene_dir, cell_type, "_Candidate_Genes.txt")
  current_candidate_peaks <- paste0(cell_type_candidate_peak_pseudo_dir, cell_type, "_Candidate_Peaks.txt")
  
  current_scRNA_cell_metadata <- paste0(cell_type_scRNA_cell_metadata_dir, cell_type, "_HVL_RNA_cell_metadata.tsv")
  current_scATAC_cell_metadata <- paste0(cell_type_scATAC_cell_metadata_dir, cell_type, "_HVL_ATAC_cell_metadata.tsv")

  current_scRNA_read_counts <- paste0(cell_type_scRNA_read_counts_dir, cell_type, "_HVL_RNA_read_counts.rds")
  current_scATAC_read_counts <- paste0(cell_type_scATAC_read_counts_dir, cell_type, "_HVL_ATAC_read_counts.rds")
  
  # Check that files exist
  # file.exists(c(hg38_ref_seq, peak_coordinates, gene_list, motifs, motif_prior, tads, current_candidate_genes, current_candidate_peaks, 
  # current_scRNA_cell_metadata, current_scATAC_cell_metadata, current_scRNA_read_counts, current_scATAC_read_counts))
  
  loaded_data <- Data_loading(current_candidate_genes, current_candidate_peaks,
                              current_scRNA_read_counts, gene_list, current_scRNA_cell_metadata,
                              current_scATAC_read_counts, peak_coordinates, current_scATAC_cell_metadata,
                              motif_prior, motifs, hg38_ref_seq)
  
  candidate_circuits <- Candidate_circuits_construction_with_TAD(loaded_data, tads)
  
  initial_model <- MAGICAL_initialization(loaded_data, candidate_circuits)
  
  circuits_linkage_posterior <- MAGICAL_estimation(loaded_data, candidate_circuits, initial_model, iteration_num = 1000)
  
  MAGICAL_circuits_output(Output_file_path = paste0(magical_output_dir, cell_type, "_MAGICAL_selected_regulatory_circuits.txt"), 
                                                    candidate_circuits, circuits_linkage_posterior)
}

# Create overall MAGICAL circuit table
cell_type <- cell_types[1]
overall_magical_df <- read.table(paste0(magical_output_dir, cell_type, "_MAGICAL_selected_regulatory_circuits.txt"), sep = "\t", header = TRUE)
overall_magical_df$Cell_Type <- cell_type
rest_of_cell_types <- cell_types[c(-1)]
for(cell_type in rest_of_cell_types) {
  current_magical_df <- read.table(paste0(magical_output_dir, cell_type, "_MAGICAL_selected_regulatory_circuits.txt"), sep = "\t", header = TRUE)
  current_magical_df$Cell_Type <- cell_type
  overall_magical_df <- rbind(overall_magical_df, current_magical_df)
}
overall_magical_df <- overall_magical_df[,c(9,1,2,3,4,5,6,7,8)]

# Create overall pseudobulk motif enrichment table
cell_type <- cell_types[1]
overall_pseudobulk_motif_enrichment_df <- read.table(paste0(sc_das_dir, "diff_peaks/", cell_type, "_D28_D1_motif_up_pseudobulk_only.tsv"), sep = "\t", header = TRUE)
overall_pseudobulk_motif_enrichment_df$Cell_Type <- cell_type
overall_pseudobulk_motif_enrichment_df$Direction <- "Pos"
current_pseudobulk_motif_enrichment_df <- read.table(paste0(sc_das_dir, "diff_peaks/", cell_type, "_D28_D1_motif_down_pseudobulk_only.tsv"), sep = "\t", header = TRUE)
current_pseudobulk_motif_enrichment_df$Cell_Type <- cell_type
current_pseudobulk_motif_enrichment_df$Direction <- "Neg"
rest_of_cell_types <- cell_types[c(-1)]
for(cell_type in rest_of_cell_types) {
  for(direction in c("up", "down")) {
    current_df <- read.table(paste0(sc_das_dir, "diff_peaks/", cell_type, "_D28_D1_motif_", direction, "_pseudobulk_only.tsv"), sep = "\t", header = TRUE)
    current_df$Cell_Type <- cell_type
    if(direction == "up") {
      current_df$Direction <- "Up"
    } else {
      current_df$Direction <- "Down"
    }
    overall_pseudobulk_motif_enrichment_df <- rbind(overall_pseudobulk_motif_enrichment_df, current_df)
  }
}
overall_pseudobulk_motif_enrichment_df <- overall_pseudobulk_motif_enrichment_df[overall_pseudobulk_motif_enrichment_df$p_adj < 0.05,]
overall_pseudobulk_motif_enrichment_df <- overall_pseudobulk_motif_enrichment_df[,c(5,6,1,2,3,4)]
overall_pseudobulk_motif_enrichment_df$TF <- sub("_.*", "", overall_pseudobulk_motif_enrichment_df$TF)

# TODO: Add N/A if not found in original sc cell types
overall_magical_df <- fill_in_info_for_magical_output(overall_magical_df, sc_das_dir, 
                                                                sc_pseudobulk_deg_combined_cell_types_table, sc_pseudobulk_deg_table,
                                                      high_pos_pseudobulk_sc_genes_bulk_passing_df$gene, 
                                                      high_neg_pseudobulk_sc_genes_bulk_passing_df$gene, sc_peaks)
write.table(overall_magical_df,
            file = paste0(magical_output_dir, "MAGICAL_overall_output.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

# For TFs associated with the current circuit, were any found to be genes in a MAGICAL circuit?
# For TFs associated with the current circuit, were any found to be DEGs in the combined cell types used as input for MAGICAL?
# For TFs associated with the current circuit, were any found to be DEGs in the more granular cell types from the original DEG analysis?
# For TFs associated with the current circuit, were any found to be DEGs in the bulk data?
# For TFs associated with the current circuit, were any found to be significant in motif enrichment analysis for the scATAC-seq data?
overall_magical_tf_df <- fill_in_info_for_magical_tf_output(overall_magical_df, overall_pseudobulk_motif_enrichment_df)
write.table(overall_magical_tf_df,
            file = paste0(magical_output_dir, "MAGICAL_overall_output_tf.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
