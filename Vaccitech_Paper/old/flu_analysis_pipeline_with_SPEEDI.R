library(SPEEDI)

################## SETUP ##################
date <- Sys.Date()
home_dir <- "~/"

# Load information about samples
sample_metadata <- read.table(paste0(home_dir, "/all_metadata_sheet.tsv"), sep = "\t", header = TRUE)

# Declare data and analysis type
data_type <- "single_cell" # Can be single_cell or multiome 
analysis_type <- "sample_paired" # Can be RNA, ATAC, sample_paired, or true_multiome

# Directory where all data are for data_type specified above
data_path <- paste0(home_dir, data_type, "/data/")
# List of samples that can potentially be processed for the data type
all_sample_id_list <- list.dirs(data_path, recursive = FALSE)
all_sample_id_list <- strsplit(all_sample_id_list, "/")
all_sample_id_list <- unlist(lapply(all_sample_id_list, tail, n = 1L))

# data_token is used to choose subset of data that we want to analyze (pre-defined in flu_data_tokens.tsv)
data_token <- "single_cell_paired_sample"
# Create directory for this particular data in the analysis directory if it doesn't exist
base_analysis_dir <- paste0(home_dir, data_type, "/analysis/", data_token, "/")
if (!dir.exists(base_analysis_dir)) {dir.create(base_analysis_dir)}
# Grab samples that we want to analyze
data_tokens <- read.table(paste0(home_dir, "flu_data_tokens.tsv"), header = TRUE)
samples <- data_tokens[data_tokens$token == data_token,]$samples
samples <-  unlist(strsplit(samples, ","))
sample_id_list <- all_sample_id_list[all_sample_id_list %in% samples]
if(length(sample_id_list) != length(samples)) {
  print(paste0("Requested samples: ", samples))
  print(paste0("Available samples (from your list): ", sample_id_list))
  stop("You have at least one invalid sample in your list")
}
sample_count <- length(sample_id_list)

# analysis_token is used to define a specific analysis
analysis_token <- "default"
analysis_dir <- paste0(home_dir, data_type, "/analysis/", data_token, "/", analysis_token, "/")
if (!dir.exists(analysis_dir)) {dir.create(analysis_dir)}
# Directory for analysis plots
plot_dir <- paste0(analysis_dir, "plots/")
if (!dir.exists(plot_dir)) {dir.create(plot_dir)}
# Directory for markers
marker_dir <- paste0(analysis_dir, "markers/")
if (!dir.exists(marker_dir)) {dir.create(marker_dir)}

# Identify high / low viral load samples, D28 / D-1 samples, and M / F samples
# Subset metadata to only contain our sample_ids
# We will add these metadata to our Seurat object later
sample_metadata <- sample_metadata[sample_metadata$aliquot_id %in% sample_id_list,]

high_viral_load_samples <- sort(sample_metadata[sample_metadata$viral_load == "high",]$aliquot)
low_viral_load_samples <- sort(sample_metadata[sample_metadata$viral_load == "low",]$aliquot)
all_viral_load_samples <- c(high_viral_load_samples, low_viral_load_samples)
d28_samples <- sort(sample_metadata[sample_metadata$time_point == "D28",]$aliquot)
d_minus_1_samples <- sort(sample_metadata[sample_metadata$time_point == "D-1",]$aliquot)
all_day_samples <- c(d28_samples, d_minus_1_samples)
male_samples <- sort(sample_metadata[sample_metadata$sex == "M",]$aliquot)
female_samples <- sort(sample_metadata[sample_metadata$sex == "F",]$aliquot)
all_sex_samples <- c(male_samples, female_samples)

# Organize sample list by high viral load followed by low viral load
sample_id_list <- sample_id_list[order(match(sample_id_list, all_viral_load_samples))]
sample_metadata <- sample_metadata[order(match(sample_metadata$aliquot,all_viral_load_samples)),]

# Sort other types of samples according to viral load order
#all_day_samples <- all_day_samples[order(match(all_day_samples, all_viral_load_samples))]
#all_sex_samples <- all_sex_samples[order(match(all_sex_samples, all_viral_load_samples))]

################## ANALYSIS ##################
if(analysis_type == "RNA" | analysis_type == "sample_paired") {
  
  # Run QC (optional)
  if(run_qc) {
    sc_obj <- generate_qc_plots(all_sc_exp_matrices, plot_dir, date, high_viral_load_samples, low_viral_load_samples,
                    d28_samples, d_minus_1_samples, male_samples, female_samples)
    # If you run QC, you probably don't want to proceed with the rest of the pipeline, so we stop here
  } else {
    # Step 2 - filter raw data (either with adaptive thresholds or strict thresholds)
    sc_obj <- FilterRawData(all_sc_exp_matrices, human = TRUE, record_doublets = FALSE, adaptive_QC_thresholds = FALSE)
    rm(all_sc_exp_matrices)
    sc_obj <- add_sample_metadata(sc_obj, high_viral_load_samples, low_viral_load_samples,
                        d28_samples, d_minus_1_samples, male_samples, female_samples)
    # Step 3 - normalize data
    sc_obj <- InitialProcessing(sc_obj, human = TRUE)
    save_seurat(save_progress, sc_obj, file = paste0(analysis_dir, "3_sc_obj.rds"))
    # Step 4 - find batches in data
    sc_obj <- InferBatches(sc_obj)
    # Step 5 - integrate data by batch
    sc_obj <- IntegrateByBatch(sc_obj)
    save_seurat(save_progress, sc_obj, file = paste0(analysis_dir, "5_sc_obj.rds"))
    # Step 6 - process integrated assay and potentially prepare to find markers via PrepSCTFindMarkers
    sc_obj <- VisualizeIntegration(sc_obj, prep_sct_find_markers = TRUE)
    save(sc_obj, file = paste0(analysis_dir, "6_sc_obj.rds"))
    # Load reference and remove cell types we don't like
    reference <- LoadReference("PBMC", human = TRUE)
    idx <- which(reference$celltype.l2 %in% c("Doublet", "B intermediate", "CD4 CTL", "gdT", "dnT", "ILC"))
    reference <- reference[,-idx]
    # Map cell types from reference to our data
    sc_obj <- MapCellTypes(sc_obj, reference, data_type = "snRNA")
    rm(reference)
    # Combine cell types and re-do majority vote
    sc_obj <- combine_cell_types_initial(sc_obj, resolution = 3)
    # Print UMAP by cell type (majority vote) and by cluster number - it will currently be messy
    print_UMAP_stage_1(sc_obj, sample_count, plot_dir, date)
    # tagged_samples <- c("3247c65ecdbfe34a", "d360f89cf9585dfe", "48ebe8475317ba95", "3c4540710e55f7b1", "fba8595c48236db8") # 14 multiome
    # tagged_samples <- c("3247c65ecdbfe34a", "d360f89cf9585dfe") # 17 multiome
    # sc_obj <- subset(x = sc_obj, subset = viral_load %in% "HVL") # HVL (10 multiome if starting with 14 multiome)
    tagged_samples <- c("48ebe8475317ba95", "3c4540710e55f7b1", "fba8595c48236db8") # 16 multiome
    # print_UMAP_tagged(sc_obj, tagged_samples, sample_count, plot_dir, date)
    # We always want to save our sc_obj after processing data through SPEEDI
    save(sc_obj, file = paste0(analysis_dir, "7_sc_obj.rds"))
    # Load sc_obj
    load(file = paste0(analysis_dir, "7_sc_obj.rds"))
    #load(paste0(analysis_dir, "singler_labels.rds"))
    # We can use clustree to help us figure out the best resolution
    # NOTE: clustree may not be that useful in the integrated setting because it'll over-cluster according to sample
    # print_clustree_plot(sc_obj, plot_dir, date)
    # Re-run majority vote with best resolution
    best_res <- 1
    sc_obj <- MajorityVote(sc_obj, best_res)
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

