library(SPEEDI)
library(Seurat)
# Load extra RNA functions
home_dir <- "~/"
source(paste0(home_dir, "extra_functions/rna/preprocessing_and_qc.R"))
source(paste0(home_dir, "extra_functions/rna/get_stats.R"))
source(paste0(home_dir, "extra_functions/rna/manipulate_data.R"))
source(paste0(home_dir, "extra_functions/rna/differential_expression.R"))
source(paste0(home_dir, "extra_functions/rna/visualization.R"))

################## SETUP ##################
date <- Sys.Date()
data_path <- "~/single_cell/data"
reference_dir <- "~/references/"
reference_file_name <- "pbmc_multimodal.h5seurat"
output_dir <- "~/single_cell/analysis/"
analysis_name <- "primary_analysis_6_subject_12_sample"
reference_tissue <- "pbmc_full"
reference_cell_type_attribute <- "celltype.l2"
species <- "human"
record_doublets <- TRUE
data_type <- "RNA"
# data_token is used to choose subset of data that we want to analyze (pre-defined in flu_data_tokens.tsv)
data_token <- "single_cell_paired_sample"
# Grab samples that we want to analyze
data_tokens <- read.table(paste0(home_dir, "flu_data_tokens.tsv"), header = TRUE)
all_sample_id_list <- list.dirs(data_path, recursive = FALSE)
all_sample_id_list <- strsplit(all_sample_id_list, "/")
all_sample_id_list <- unlist(lapply(all_sample_id_list, tail, n = 1L))
samples <- data_tokens[data_tokens$token == data_token,]$samples
samples <-  unlist(strsplit(samples, ","))
sample_id_list <- all_sample_id_list[all_sample_id_list %in% samples]

# Load information about samples
sample_metadata <- read.table(paste0(home_dir, "/all_metadata_sheet.tsv"), sep = "\t", header = TRUE)
sample_metadata <- sample_metadata[sample_metadata$aliquot_id %in% sample_id_list,]
sample_metadata_for_SPEEDI_df <- sample_metadata
rownames(sample_metadata_for_SPEEDI_df) <- sample_metadata_for_SPEEDI_df$aliquot_id
sample_metadata_for_SPEEDI_df <- sample_metadata_for_SPEEDI_df[c("subject_id", "time_point", "sex", "viral_load")]
sample_metadata_for_SPEEDI_df$time_point[sample_metadata_for_SPEEDI_df$time_point == 'D-1'] <- 'D_minus_1'

# Normalize paths (in case user provides relative paths)
data_path <- normalize_dir_path(data_path)
reference_dir <- normalize_dir_path(reference_dir)
output_dir <- normalize_dir_path(output_dir)
# ArchR likes to write some files to the working directory, so we'll set our working directory to output_dir
# and then reset it to the original working directory once we're done running SPEEDI
old_wd <- getwd()
# Create output_dir if it doesn't already exist
if (!dir.exists(output_dir)) {dir.create(output_dir)}
# Add "/" to end of output_dir if not already present
last_char_of_output_dir_path <- substr(output_dir, nchar(output_dir), nchar(output_dir))
if(last_char_of_output_dir_path != "/") {
  output_dir <- paste0(output_dir, "/")
}
# Set analysis name
if(is.null(analysis_name)) {
  analysis_name <- paste0(gsub(" ", "_", Sys.time()), "_SPEEDI")
  analysis_name <- gsub(":", "-", analysis_name)
}
# Update our output dir to be the specific analysis directory
output_dir <- paste0(output_dir, analysis_name, "/")
if (!dir.exists(output_dir)) {dir.create(output_dir)}
setwd(output_dir)
# Create log file
log_file_name <- paste0(gsub(" ", "_", Sys.time()), "_SPEEDI")
log_file_name <- gsub(":", "-", log_file_name)
log_file_name <- paste0(output_dir, log_file_name)
log_file <- logr::log_open(log_file_name, logdir = FALSE)
# Load reference
reference <- LoadReferenceSPEEDI(reference_tissue = reference_tissue, species = species, reference_dir = reference_dir,
                                 reference_file_name = reference_file_name, log_flag = TRUE)
idx <- which(reference$celltype.l2 %in% c("Doublet", "B intermediate", "CD4 CTL", "gdT", "dnT", "ILC"))
reference <- reference[,-idx]

# Output dirs for RNA and ATAC
RNA_output_dir <- paste0(output_dir, "RNA", "/")
if (!dir.exists(RNA_output_dir)) {dir.create(RNA_output_dir)}
# Read in RNA data, filter data, perform initial processing, infer batches, integrate by batch, and process UMAP of integration
all_sc_exp_matrices <- Read_RNA(data_path = data_path, sample_id_list = sample_id_list, log_flag = TRUE)
sc_obj <- FilterRawData_RNA(all_sc_exp_matrices = all_sc_exp_matrices, species = species,
                            record_doublets = record_doublets, output_dir = RNA_output_dir,
                            log_file_path = log_file_name, log_flag = TRUE)
rm(all_sc_exp_matrices)
sc_obj <- InitialProcessing_RNA(sc_obj = sc_obj, species = species, metadata_df = sample_metadata_for_SPEEDI_df, log_flag = TRUE)
sc_obj <- InferBatches_alt(sc_obj = sc_obj, log_flag = TRUE) # STOPPED AFTER THIS STEP - 6 batches with new approach instead of 3? Weird? 
# save(sc_obj, file = paste0(RNA_output_dir, analysis_name, ".RNA.rds"))
# load(paste0(RNA_output_dir, "primary_analysis_6_subject_12_sample.RNA_old.rds"))
sc_obj <- IntegrateByBatch_RNA(sc_obj = sc_obj, log_flag = TRUE)
sc_obj <- VisualizeIntegration(sc_obj = sc_obj, log_flag = TRUE)
sc_obj <- MapCellTypes_RNA(sc_obj = sc_obj, reference = reference,
                           reference_cell_type_attribute = reference_cell_type_attribute,
                           output_dir = RNA_output_dir, log_flag = TRUE)
# sc_obj <- MajorityVote_RNA_alt(sc_obj)
save(sc_obj, file = paste0(RNA_output_dir, analysis_name, ".RNA.old.algorithm.rds"))
# load(paste0(RNA_output_dir, "primary_analysis_6_subject_12_sample.RNA.old.algorithm.rds"))

sc_obj$old.predicted.id <- sc_obj$predicted.id
Cell_type_combined <- sc_obj$predicted.id
idx <- grep("CD4 T", Cell_type_combined)
Cell_type_combined[idx] <- "CD4 Memory"
idx <- grep("CD8 T", Cell_type_combined)
Cell_type_combined[idx] <- "CD8 Memory"
idx <- grep("cDC", Cell_type_combined)
Cell_type_combined[idx] <- "cDC"
idx <- grep("Proliferating", Cell_type_combined)
Cell_type_combined[idx] <- "Proliferating"
idx <- grep("CD4 Naive", Cell_type_combined)
Cell_type_combined[idx] <- "T Naive"
idx <- grep("CD8 Naive", Cell_type_combined)
Cell_type_combined[idx] <- "T Naive"
idx <- grep("Treg", Cell_type_combined)
Cell_type_combined[idx] <- "T Naive"
sc_obj$predicted.id <- Cell_type_combined
sc_obj <- MajorityVote_RNA_alt(sc_obj)

print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Majority_Vote_Cell_Type.png",
               group_by_category = "predicted_celltype_majority_vote", output_dir = output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Cluster.png",
               group_by_category = "seurat_clusters", output_dir = output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Raw_Predicted_Cell_Type.png",
               group_by_category = "predicted.id", output_dir = output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Sample.png",
               group_by_category = "sample", output_dir = output_dir,
               log_flag = log_flag)







combined_cell_type_dir <- paste0(RNA_output_dir, "combined_cell_types/")

print_UMAP_RNA(sc_obj, paste0("pre.clusters_by_cell_type_majority_vote_", date, ".png"), "predicted_celltype_majority_vote", combined_cell_type_dir)
print_UMAP_RNA(sc_obj, paste0("pre.clusters_by_cluster_num_", date, ".png"), "seurat_clusters", combined_cell_type_dir)
print_UMAP_RNA(sc_obj, paste0("pre.clusters_by_cell_type_", date, ".png"), "predicted.id", combined_cell_type_dir)
print_UMAP_RNA(sc_obj, paste0("pre.clusters_by_viral_load_", date, ".png"), "viral_load", combined_cell_type_dir)
print_UMAP_RNA(sc_obj, paste0("pre.clusters_by_sample_", date, ".png"), "sample", combined_cell_type_dir)










    # sc_obj <- subset(x = sc_obj, subset = viral_load %in% "HVL") # HVL (10 multiome if starting with 14 multiome)
    # We always want to save our sc_obj after processing data through SPEEDI
    save(sc_obj, file = paste0(analysis_dir, "7_sc_obj.rds"))
    # Load sc_obj
    load(file = paste0(analysis_dir, "7_sc_obj.rds"))
    #load(paste0(analysis_dir, "singler_labels.rds"))
    
    # To decide which clusters we need to remove, we will capture information about clusters
    # We will also run DE for each cluster to find cell type markers
    #raw_cluster_info <- capture_cluster_info(sc_obj)
    #run_differential_expression_cluster(sc_obj, marker_dir)
    # Remove messy clusters
    messy_clusters <- c(6,7,13,14,15,19,28) # 14 sample multiome
    #messy_clusters <- c(2,13,16,18,19,20,21,24,28,32) # 19 sample multiome
    idxPass <- which(Idents(sc_obj) %in% messy_clusters)
    cellsPass <- names(sc_obj$orig.ident[-idxPass])
    sc_obj.minus.messy.clusters <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
    # Override labels manually where necessary
    #sc_obj.minus.messy.clusters <- override_cluster_label(sc_obj.minus.messy.clusters, c(15), "CD16 Mono")
    # sc_obj.minus.messy.clusters <- override_cluster_label(sc_obj.minus.messy.clusters, c(25), "pDC") # 19 sample multiome
    sc_obj.minus.messy.clusters <- override_cluster_label(sc_obj.minus.messy.clusters, c(20), "CD4 Memory") # 14 sample multiome?
    # Remove specific samples that we're not interested in
    # sc_obj.minus.messy.clusters <- remove_specific_samples_from_sc_obj(sc_obj.minus.messy.clusters, tagged_samples)
    # print_UMAP(sc_obj.minus.messy.clusters.removed.samples, sample_count, "tagged", plot_dir, paste0("post.clusters_all_untagged_", date, ".png"))
    # Print info about sample representation and breakdown of categories per cell type
    #print(table(sc_obj.minus.messy.clusters$sample))
    #print_celltype_counts(sc_obj.minus.messy.clusters)
    # Remove noisy cells
    sc_obj.minus.messy.clusters <- remove_cells_based_on_umap(sc_obj.minus.messy.clusters, -2.25, 0, 1, 2.5) # 14 sample multiome
    sc_obj.minus.messy.clusters <- remove_cells_based_on_umap(sc_obj.minus.messy.clusters, 1.5, 2.5, -3, -1.5) # 14 sample multiome
    sc_obj.minus.messy.clusters <- remove_cells_based_on_umap(sc_obj.minus.messy.clusters, 2.2, 3.5, -6, -4) # 14 sample multiome
    #sc_obj.minus.messy.clusters <- remove_cells_based_on_umap(sc_obj.minus.messy.clusters, -4, -2, 1, 6) # 19 sample multiome
    #sc_obj.minus.messy.clusters <- remove_cells_based_on_umap(sc_obj.minus.messy.clusters, 1, 3, -2, 0) # 19 sample multiome
    #sc_obj.minus.messy.clusters <- remove_cells_based_on_umap(sc_obj.minus.messy.clusters, 3, 4, -3, -1.5) # 19 sample multiome
    # Print plots for various metadata categories for final UMAP
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "predicted_celltype_majority_vote", plot_dir, paste0("post.clusters_by_cell_type_majority_vote_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "seurat_clusters", plot_dir, paste0("post.clusters_by_cluster_num_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "predicted.id", plot_dir, paste0("post.clusters_by_cell_type_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "viral_load", plot_dir, paste0("post.clusters_by_viral_load_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "sample", plot_dir, paste0("post.clusters_by_sample_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "day", plot_dir, paste0("post.clusters_by_day_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "sex", plot_dir, paste0("post.clusters_by_sex_", date, ".png"))
    # write cells and associated cell type majority predictions to file (for ATAC-seq labeling)
    cells_for_ATAC <- data.frame("cells" = sc_obj.minus.messy.clusters$cell_name, voted_type = sc_obj.minus.messy.clusters$predicted_celltype_majority_vote)
    write.csv(cells_for_ATAC, file = paste0(analysis_dir, "rna_seq_labeled_cells_", date, "-14_final.csv"), quote = FALSE, row.names = FALSE)
    # Combine cell types for MAGICAL and other analyses that require ATAC-seq (granularity isn't as good for ATAC-seq)
    sc_obj.minus.messy.clusters <- combine_cell_types_magical(sc_obj.minus.messy.clusters)
    # Run differential expression for each cell type within each group of interest
    differential_genes_dir <- paste0(analysis_dir, "diff_genes/", date, "/final/")
    if (!dir.exists(differential_genes_dir)) {dir.create(differential_genes_dir, recursive = TRUE)}
    run_differential_expression_group(sc_obj.minus.messy.clusters, differential_genes_dir, "viral_load")
    run_differential_expression_group(sc_obj.minus.messy.clusters, differential_genes_dir, "day")
    run_differential_expression_group(sc_obj.minus.messy.clusters, differential_genes_dir, "sex")
    create_magical_cell_type_proportion_file(sc_obj.minus.messy.clusters, "viral_load", high_viral_load_samples, d28_samples, male_samples)
    create_magical_cell_type_proportion_file(sc_obj.minus.messy.clusters, "day", high_viral_load_samples, d28_samples, male_samples, token = "14_final")
    create_magical_cell_type_proportion_file(sc_obj.minus.messy.clusters, "sex", high_viral_load_samples, d28_samples, male_samples)
    pseudobulk_rna_dir <- paste0(analysis_dir, "pseudobulk_rna/", date, "/")
    if (!dir.exists(pseudobulk_rna_dir)) {dir.create(pseudobulk_rna_dir, recursive = TRUE)}
    create_magical_cell_type_pseudobulk_files(sc_obj.minus.messy.clusters, pseudobulk_rna_dir, token = "14_final")
  }
} else if(analysis_type == "ATAC_seq") {
  # Label input files
  inputFiles <- paste0(data_path, sample_id_list, "/outs/atac_fragments.tsv.gz")
  names(inputFiles) <- all_viral_load_samples
  # Add relevant genome for ArchR
  addArchRGenome("hg38")
  # Create ArchR project from input files
  proj <- load_archR_from_input_files(inputFiles, analysis_dir)
  # Re-load ArchR project
  proj <- loadArchRProject(path = paste0(analysis_dir, "/ArchR/"))
  # Add sample metadata
  proj <- add_sample_metadata_atac(proj, high_viral_load_samples, low_viral_load_samples,
                                   d28_samples, d_minus_1_samples, male_samples, female_samples)
  viral_load_metadata <- parse_metadata_for_samples(proj, "viral_load", high_viral_load_samples, low_viral_load_samples,
                                                    d28_samples, d_minus_1_samples, male_samples, female_samples)
  day_metadata <- parse_metadata_for_samples(proj, "day", high_viral_load_samples, low_viral_load_samples,
                                             d28_samples, d_minus_1_samples, male_samples, female_samples)
  sex_metadata <- parse_metadata_for_samples(proj, "sex", high_viral_load_samples, low_viral_load_samples,
                                             d28_samples, d_minus_1_samples, male_samples, female_samples)
  # Plot what dataset looks like before any processing
  plot_qc_atac(proj, date)
  # Filter out cells that don't meet TSS enrichment / doublet enrichment / nucleosome ratio criteria
  idxPass <- which(proj$TSSEnrichment >= 8 & proj$NucleosomeRatio < 2 & proj$DoubletEnrichment < 3)
  cellsPass <- proj$cellNames[idxPass]
  proj <- proj[cellsPass, ]
  # List number of cells remaining for each condition and each sample
  table(proj$viral_load)
  table(proj$day)
  table(proj$sex)
  table(proj$Sample)
  # Perform dimensionality reduction on cells (addIterativeLSI), create UMAP embedding (addUMAP),
  # and add cluster information (addClusters)
  proj <- dimensionality_reduc(proj)
  plot_atac_after_filtering(proj, date)
  # Map from scRNA reference to ATAC data
  scRNA_reference <- load_rna_reference_for_atac("~/reference/")
  proj <- map_reference_to_atac(proj)
  rm(scRNA_reference)
  plot_atac_after_integration(proj, date)
  saveArchRProject(ArchRProj = proj, load = FALSE)
  # Load ArchR project
  proj <- loadArchRProject(path = paste0(analysis_dir, "/ArchR/"))
  if(use_rna_labels) {
    proj <- add_rna_labels_for_atac_data(proj, analysis_dir, source_rna_file = "rna_seq_labeled_cells_2023-04-10-14_final.csv", subset_to_rna)
  }
  proj <- combine_cell_types_atac(proj)
  # If we subset to RNA, we don't need to do any majority voting in clusters, etc.
  # Otherwise, we do!
  if(subset_to_rna) {
    final_proj <- proj
    # Just for convenience in code below
    final_proj$Cell_type_voting <- final_proj$predictedGroup
    final_proj <- remove_cell_types(final_proj, c("HSPC", "Plasmablast", "Proliferating", "Platelet", "MAIT"))
    final_proj <- remove_cells_based_on_umap_atac(final_proj, -2, 1, -3, 2) # Multiome 14
    final_proj <- remove_cells_based_on_umap_atac(final_proj, 1, 4.5, -4.5, -2.5) # Multiome 14
    final_proj <- remove_cells_based_on_umap_atac(final_proj, -2.5, -2, -2.5, 0.5) # Multiome 14
    #final_proj <- remove_cells_based_on_umap_atac(final_proj, -3.5, 1, 0, 3.5) # Multiome 19
    #final_proj <- remove_cells_based_on_umap_atac(final_proj, 1, 5, 2.75, 5) # Multiome 19
    plot_atac_after_majority_vote_or_subset(final_proj, date)
  } else {
    proj <- perform_majority_vote(proj)
    plot_atac_after_majority_vote_or_subset(proj, date)
    proj <- remove_single_cell_clusters(proj)
    cluster_info <- get_cluster_info(proj)
    proj <- override_cluster_label_atac(proj, c("C8"), "CD16 Mono")
    #Remove the messy clusters (determined through visual inspection and seeing distribution of cells in each cluster)
    idxPass <- which(proj$Clusters %in% c("C1", "C29", "C30", "C31", "C38"))
    cellsPass <- proj$cellNames[-idxPass]
    final_proj <- proj[cellsPass, ]
    plot_atac_after_majority_vote_or_subset(final_proj, date)
  }
  print_cell_type_distributions(final_proj)
  create_cell_type_proportion_MAGICAL_atac(final_proj, analysis_dir, c("day"), day_metadata)
  addArchRGenome("hg38")
  final_proj <- pseudo_bulk_replicates_and_call_peaks(final_proj)
  # Create peak matrix (matrix containing insertion counts within our merged peak set) for differential accessibility
  # calculations
  final_proj <- addPeakMatrix(final_proj)
  differential_peaks_dir <- paste0(analysis_dir, "diff_peaks-14_final/", date, "/")
  if (!dir.exists(differential_peaks_dir)) {dir.create(differential_peaks_dir, recursive = TRUE)}
  calculate_daps_for_each_cell_type(final_proj, differential_peaks_dir)
  # Create Peaks.txt file for MAGICAL
  peak_txt_file <- create_peaks_file(final_proj, analysis_dir)
  # Create peak_motif_matches.txt file for MAGICAL
  create_peak_motif_matches_file(final_proj, analysis_dir, peak_txt_file)
  # Create pseudobulk counts for peaks for each cell type
  pseudo_bulk_dir <- paste0(analysis_dir, "pseudo_bulk_atac-14_final/", date, "/")
  if (!dir.exists(pseudo_bulk_dir)) {dir.create(pseudo_bulk_dir, recursive = TRUE)}
  create_pseudobulk_atac(final_proj, pseudo_bulk_dir)
} else {
  stop("Invalid analysis type")
}

