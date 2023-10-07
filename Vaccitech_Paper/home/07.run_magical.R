# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

source(paste0(base_dir, "extra_functions/MAGICAL_functions.R"))

# STEP 0: Set up directories and global files
cell_type_candidate_gene_dir <- paste0(sc_magical_dir, "Candidate_Genes/")
if (!dir.exists(cell_type_candidate_gene_dir)) {dir.create(cell_type_candidate_gene_dir, recursive = TRUE)}
cell_type_candidate_peak_dir <- paste0(sc_magical_dir, "Candidate_Peaks/")
if (!dir.exists(cell_type_candidate_peak_dir)) {dir.create(cell_type_candidate_peak_dir, recursive = TRUE)}
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

# b) Cell type candidate peaks.txt

for(cell_type in unique(sc_das_lenient$Cell_Type)) {
  cell_type_sc_das_lenient <- sc_das_lenient[sc_das_lenient$Cell_Type == cell_type,]
  cell_type_sc_das_lenient <- cell_type_sc_das_lenient[,c("chr", "start", "end")]
  write.table(cell_type_sc_das_lenient,
              file = paste0(cell_type_candidate_peak_dir, sub(" ", "_", cell_type), "_Candidate_Peaks.txt"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

# STEP 2: FIND MATCHING CELL TYPES FROM RNA / ATAC FOR MAGICAL ANALYSIS
possible_RNA_types <- list.files(path = cell_type_candidate_gene_dir)
possible_RNA_types <- unlist(strsplit(possible_RNA_types, "_Candidate_Genes.txt"))

possible_ATAC_types <- list.files(path = cell_type_candidate_peak_dir)
possible_ATAC_types <- unlist(strsplit(possible_ATAC_types, "_Candidate_Peaks.txt"))

cell_types <- intersect(possible_RNA_types, possible_ATAC_types)

# STEP 3: RUN MAGICAL
set.seed(get_speedi_seed())
for(cell_type in cell_types) {
  print(cell_type)
  current_candidate_genes <- paste0(cell_type_candidate_gene_dir, cell_type, "_Candidate_Genes.txt")
  current_candidate_peaks <- paste0(cell_type_candidate_peak_dir, cell_type, "_Candidate_Peaks.txt")
  
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
  
  MAGICAL_circuits_output(Output_file_path = paste0(magical_output_dir, cell_type, "_MAGICAL_selected_regulatory_circuits_6.txt"), 
                                                    candidate_circuits, circuits_linkage_posterior)
}
