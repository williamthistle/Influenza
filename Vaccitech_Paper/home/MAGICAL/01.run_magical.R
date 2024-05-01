### THIS FUNCTION:
### RUNS MAGICAL AND CREATES 1) OVERALL MAGICAL DF, AND 2) TF-CENTRIC MAGICAL DF

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

source(paste0(base_dir, "extra_functions/MAGICAL_functions.R"))

# STEP 0: Set up directories and global files
cell_type_candidate_gene_dir <- paste0(MAGICAL_hvl_placebo_dir, "Candidate_Genes/")
if (!dir.exists(cell_type_candidate_gene_dir)) {dir.create(cell_type_candidate_gene_dir, recursive = TRUE)}
cell_type_candidate_peak_dir <- paste0(MAGICAL_hvl_placebo_dir, "Candidate_Peaks/")
if (!dir.exists(cell_type_candidate_peak_dir)) {dir.create(cell_type_candidate_peak_dir, recursive = TRUE)}
cell_type_scATAC_cell_metadata_dir <- paste0(MAGICAL_hvl_placebo_dir, "scATAC_Cell_Metadata/")
if (!dir.exists(cell_type_scATAC_cell_metadata_dir)) {dir.create(cell_type_scATAC_cell_metadata_dir, recursive = TRUE)}
cell_type_scRNA_cell_metadata_dir <- paste0(MAGICAL_hvl_placebo_dir, "scRNA_Cell_Metadata/")
if (!dir.exists(cell_type_scRNA_cell_metadata_dir)) {dir.create(cell_type_scRNA_cell_metadata_dir, recursive = TRUE)}
cell_type_scATAC_read_counts_dir <- paste0(MAGICAL_hvl_placebo_dir, "scATAC_Read_Counts/")
if (!dir.exists(cell_type_scATAC_read_counts_dir)) {dir.create(cell_type_scATAC_read_counts_dir, recursive = TRUE)}
cell_type_scRNA_read_counts_dir <- paste0(MAGICAL_hvl_placebo_dir, "scRNA_Read_Counts/")
if (!dir.exists(cell_type_scRNA_read_counts_dir)) {dir.create(cell_type_scRNA_read_counts_dir, recursive = TRUE)}
magical_output_dir <- paste0(MAGICAL_hvl_placebo_dir, "Output/")
if (!dir.exists(magical_output_dir)) {dir.create(magical_output_dir, recursive = TRUE)}

hg38_ref_seq <- paste0(MAGICAL_hvl_placebo_dir, "hg38_Refseq.txt")
peak_coordinates <- paste0(MAGICAL_hvl_placebo_dir, "scATAC_Peak_Coordinates/HVL_ATAC_peak_coordinates.tsv")
gene_list <- paste0(MAGICAL_hvl_placebo_dir, "scRNA_Genes/HVL_RNA_genes.tsv")
motifs <- paste0(MAGICAL_hvl_placebo_dir, "Motifs.txt")
motif_prior <- paste0(MAGICAL_hvl_placebo_dir, "scATAC_Motif_Mapping_Prior/HVL_ATAC_motif_mapping_prior.tsv")
tads <- paste0(MAGICAL_hvl_placebo_dir, "RaoGM12878_40kb_TopDomTADs_filtered_hg38.txt")
distance_control=5e5

cell_types <- c("B", "CD14 Mono", "CD16 Mono", "Proliferating", "NK", "CD4 Memory", "CD8 Memory", "MAIT", "CD4 Naive", "CD8 Naive", "cDC")

# STEP 1: CREATE CANDIDATE GENE AND PEAK FILES

current_candidate_deg_table <- scRNA_hvl_placebo_MAGICAL_cell_types_degs
current_candidate_deg_table <- current_candidate_deg_table[current_candidate_deg_table$Cell_Type != "T Naive",]
scRNA_hvl_placebo_degs_naive_t_cells <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == "CD4 Naive" | scRNA_hvl_placebo_degs$Cell_Type == "CD8 Naive",]
current_candidate_deg_table <- rbind(current_candidate_deg_table, scRNA_hvl_placebo_degs_naive_t_cells)

# a) Cell type candidate genes.txt

for(cell_type in unique(current_candidate_deg_table$Cell_Type)) {
  write.table(current_candidate_deg_table[current_candidate_deg_table$Cell_Type == cell_type,]$Gene_Name,
              file = paste0(cell_type_candidate_gene_dir, sub(" ", "_", cell_type), "_Candidate_Genes.txt"), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}

# b) Cell type candidate peaks.txt
for(cell_type in cell_types) {
  current_das <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", sub(" ", "_", cell_type), "-time_point-controlling_for_subject_id_sc_pct_0.01.tsv"), sep = "\t",
                                       header = TRUE)
  peak_list <- rownames(current_das)
  chromosomes <- sapply(strsplit(peak_list, "-"), `[`, 1)
  start_coords <- sapply(strsplit(peak_list, "-"), `[`, 2)
  end_coords <- sapply(strsplit(peak_list, "-"), `[`, 3)
  final_das <- data.frame(chr = chromosomes, start = start_coords, end = end_coords)
  write.table(final_das,
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

# Read MAGICAL df
# overall_magical_df <- read.table(paste0(MAGICAL_hvl_placebo_output_dir, "MAGICAL_overall_output.tsv"), sep = "\t", header = TRUE)

# Add metadata - stage 1
overall_magical_df <- add_info_stage_1_MAGICAL(overall_magical_df, current_candidate_deg_table)
write.table(overall_magical_df, file = paste0(MAGICAL_hvl_placebo_output_dir, "MAGICAL_overall_output_with_stage_1_info.tsv"),
            sep = "\t", quote = FALSE)

overall_magical_df_with_pseudobulk_correction <- overall_magical_df[overall_magical_df$sc_p_val < 0.05,]
overall_magical_df_with_pseudobulk_correction <- overall_magical_df_with_pseudobulk_correction[overall_magical_df_with_pseudobulk_correction$pseudobulk_p_val < 0.05 | overall_magical_df_with_pseudobulk_correction$pseudobulk_robust_p_val < 0.05,]

# Two new tables!
magical_gene_overlap_df <- create_magical_gene_overlap_df(overall_magical_df, hvl_placebo_LRT_analysis_results_filtered)
magical_site_overlap_df <- create_magical_site_overlap_df(overall_magical_df)

# Find overlapping circuits between cell types
overlapping_circuits_df <- find_overlapping_circuits(overall_magical_df)
overlapping_circuit_cell_type_matrix <- find_overlapping_circuit_cell_type_matrix(overlapping_circuits_df)

# Convert matrix to data frame
overlapping_circuit_cell_type_df <- melt(overlapping_circuit_cell_type_matrix)

# Create the heatmap
ggplot(overlapping_circuit_cell_type_df, aes(x = Var2, y = Var1, fill = value, label = value)) +
  geom_tile(color = "white") +
  geom_text(size = 3) +  # Add text labels
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# TODO: Create UpsetPlot for genes like Aliza paper
# 

# Create targets for each cell type
magical_tf_vs_cell_type_df <- create_tf_vs_cell_type_df(overall_magical_df)
# Are these TFs found in scRNA data?
tfs_found_in_scRNA_df <- find_tfs_in_scRNA_data(magical_tf_vs_cell_type_df, current_candidate_deg_table)

magical_tf_heatmap_plot <- create_tf_heatmap_plot(tfs_found_in_scRNA_df)

# TF analysis for min.pct = 0.05
overall_magical_df_0.05 <- overall_magical_df[overall_magical_df$pct.1 >= 0.05 | overall_magical_df$pct.2 >= 0.05,]
magical_tf_vs_cell_type_df_0.05 <- create_tf_vs_cell_type_df(overall_magical_df_0.05)
tfs_found_in_scRNA_df_0.05 <- find_tfs_in_scRNA_data(magical_tf_vs_cell_type_df_0.05, current_candidate_deg_table)

# TF analysis for min.pct = 0.1
overall_magical_df_0.1 <- overall_magical_df[overall_magical_df$pct.1 >= 0.1 | overall_magical_df$pct.2 >= 0.1,]
magical_tf_vs_cell_type_df_0.1 <- create_tf_vs_cell_type_df(overall_magical_df_0.1)
tfs_found_in_scRNA_df_0.1 <- find_tfs_in_scRNA_data(magical_tf_vs_cell_type_df_0.1, current_candidate_deg_table)

# Create TF heatmap plot
magical_tf_heatmap_plot <- create_tf_heatmap_plot(tfs_found_in_scRNA_df_0.1)


# Create targets for each cell type
magical_tf_vs_cell_type_df <- create_tf_vs_cell_type_df(overall_magical_df_with_pseudobulk_correction)
# Are these TFs found in scRNA data?
tfs_found_in_scRNA_df <- find_tfs_in_scRNA_data(magical_tf_vs_cell_type_df, current_candidate_deg_table)

# TF analysis for min.pct = 0.05
overall_magical_df_with_pseudobulk_correction_0.05 <- overall_magical_df_with_pseudobulk_correction[overall_magical_df_with_pseudobulk_correction$pct.1 >= 0.05 | overall_magical_df_with_pseudobulk_correction$pct.2 >= 0.05,]
magical_tf_vs_cell_type_df_0.05 <- create_tf_vs_cell_type_df(overall_magical_df_with_pseudobulk_correction_0.05)
tfs_found_in_scRNA_df_0.05 <- find_tfs_in_scRNA_data(magical_tf_vs_cell_type_df_0.05, current_candidate_deg_table)

# TF analysis for min.pct = 0.1
overall_magical_df_with_pseudobulk_correction_0.1 <- overall_magical_df_with_pseudobulk_correction[overall_magical_df_with_pseudobulk_correction$pct.1 >= 0.1 | overall_magical_df_with_pseudobulk_correction$pct.2 >= 0.1,]
magical_tf_vs_cell_type_df_0.1 <- create_tf_vs_cell_type_df(overall_magical_df_with_pseudobulk_correction_0.1)
tfs_found_in_scRNA_df_0.1 <- find_tfs_in_scRNA_data(magical_tf_vs_cell_type_df_0.1, current_candidate_deg_table)

# Create TF heatmap plot
magical_tf_heatmap_plot <- create_tf_heatmap_plot(tfs_found_in_scRNA_df_0.1)




# Create overall motif enrichment table
fc_0.1_motif_cell_types <- c("B")
fc_0.585_motif_cell_types <- c("CD4 Memory", "CD8 Memory")
fc_1_motif_cell_types <- c("CD14 Mono", "CD16 Mono", "T Naive", "NK")
# not parsing MAIT or Proliferating because everything is significant 

cell_type <- fc_0.1_motif_cell_types[1]
overall_motif_enrichment_df <- read_motif_table(current_dir = paste0(sc_das_dir, "motif_enrichment/", sub(" ", "_", cell_type), "/sc_filtered/0.01/with_bg/"), analysis_type = "sc_filtered", pct = 0.01, fc_threshold = 0.1, direction = "pos")
overall_motif_enrichment_df$Cell_Type <- cell_type
overall_motif_enrichment_df$Direction <- "Up"
overall_motif_enrichment_df$rank <- seq.int(nrow(overall_motif_enrichment_df))
current_motif_enrichment_df <- read_motif_table(current_dir = paste0(sc_das_dir, "motif_enrichment/", sub(" ", "_", cell_type), "/sc_filtered/0.01/with_bg/"), analysis_type = "sc_filtered", pct = 0.01, fc_threshold = 0.1, direction = "neg")
current_motif_enrichment_df$Cell_Type <- cell_type
current_motif_enrichment_df$Direction <- "Down"
current_motif_enrichment_df$rank <- seq.int(nrow(current_motif_enrichment_df))
overall_motif_enrichment_df <- rbind(overall_motif_enrichment_df, current_motif_enrichment_df)


for(cell_type in fc_0.585_motif_cell_types) {
  for(direction in c("pos", "neg")) {
    current_motif_enrichment_df <- read_motif_table(current_dir = paste0(sc_das_dir, "motif_enrichment/", sub(" ", "_", cell_type), "/sc_filtered/0.01/with_bg/"), analysis_type = "sc_filtered", pct = 0.01, fc_threshold = 0.585, direction = direction)
    current_motif_enrichment_df$Cell_Type <- cell_type
    if(direction == "pos") {
      current_motif_enrichment_df$Direction <- "Up"
    } else {
      current_motif_enrichment_df$Direction <- "Down"
    }
    current_motif_enrichment_df$rank <- seq.int(nrow(current_motif_enrichment_df))
    overall_motif_enrichment_df <- rbind(overall_motif_enrichment_df, current_motif_enrichment_df)
  }
}

for(cell_type in fc_1_motif_cell_types) {
  for(direction in c("pos", "neg")) {
    current_motif_enrichment_df <- read_motif_table(current_dir = paste0(sc_das_dir, "motif_enrichment/", sub(" ", "_", cell_type), "/sc_filtered/0.01/with_bg/"), analysis_type = "sc_filtered", pct = 0.01, fc_threshold = 1, direction = direction)
    current_motif_enrichment_df$Cell_Type <- cell_type
    if(direction == "pos") {
      current_motif_enrichment_df$Direction <- "Up"
    } else {
      current_motif_enrichment_df$Direction <- "Down"
    }
    current_motif_enrichment_df$rank <- seq.int(nrow(current_motif_enrichment_df))
    overall_motif_enrichment_df <- rbind(overall_motif_enrichment_df, current_motif_enrichment_df)
  }
}

overall_motif_enrichment_df <- overall_motif_enrichment_df[overall_motif_enrichment_df$p.adjust < 0.05,]

overall_motif_enrichment_df <- overall_motif_enrichment_df[,c(10,11,8,6,7,9,12)]

overall_magical_df <- fill_in_info_for_magical_output(overall_magical_df, sc_das_dir, 
                                                                sc_pseudobulk_deg_combined_cell_types_table, sc_pseudobulk_deg_table,
                                                      high_pos_pseudobulk_sc_genes_bulk_passing_df$gene, 
                                                      high_neg_pseudobulk_sc_genes_bulk_passing_df$gene, sc_peaks, overall_motif_enrichment_df)
write.table(overall_magical_df,
            file = paste0(magical_output_dir, "MAGICAL_overall_output.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

# For TFs associated with the current circuit, were any found to be genes in a MAGICAL circuit?
# For TFs associated with the current circuit, were any found to be DEGs in the combined cell types used as input for MAGICAL?
# For TFs associated with the current circuit, were any found to be DEGs in the more granular cell types from the original DEG analysis?
# For TFs associated with the current circuit, were any found to be DEGs in the bulk data?
# For TFs associated with the current circuit, were any found to be significant in motif enrichment analysis for the scATAC-seq data?
overall_magical_tf_df <- fill_in_info_for_magical_tf_output(overall_magical_df, overall_motif_enrichment_df)
write.table(overall_magical_tf_df,
            file = paste0(magical_output_dir, "MAGICAL_overall_output_tf.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
