library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(clustree)
library(dplyr)
library(ggplot2)
library(hexbin)
library(mclust)
library(openxlsx)
library(outliers)
library(pheatmap)
library(scDblFinder)
library(SeuratDisk)
library(BiocParallel)
library(stringr)
library(writexl)
library(SoupX)
library(SingleR)

################## SETUP ##################
date <- Sys.Date()
home_dir <- "~/"
# Load SPEEDI (for RNA-seq analyses)
SPEEDI_dir <- paste0(home_dir, "SPEEDI")
source(paste0(SPEEDI_dir, "/rna/prototype_API.R"))
source(paste0(SPEEDI_dir, "/rna/preprocessing_and_qc.R"))
source(paste0(SPEEDI_dir, "/rna/get_stats.R"))
source(paste0(SPEEDI_dir, "/rna/manipulate_data.R"))
source(paste0(SPEEDI_dir, "/rna/differential_expression.R"))
source(paste0(SPEEDI_dir, "/rna/visualization.R"))
source(paste0(SPEEDI_dir, "/atac/preprocessing_and_qc.R"))
source(paste0(SPEEDI_dir, "/atac/processing.R"))

# Load information about samples
sample_metadata <- read.table(paste0(SPEEDI_dir, "/sample_metadata.tsv"), sep = "\t", header = TRUE)

# Declare data and analysis type
data_type <- "multiome" # Can be multiome or single_cell
analysis_type <- "RNA-seq" # Can be RNA-seq or ATAC-seq

# Directory where all data are for data_type specified above
data_path <- paste0(home_dir, data_type, "/data/")
# List of samples that can potentially be processed for the data type
all_sample_id_list <- list.dirs(data_path, recursive = FALSE)
all_sample_id_list <- strsplit(all_sample_id_list, "/")
all_sample_id_list <- unlist(lapply(all_sample_id_list, tail, n = 1L))

# data_token is used to choose subset of data that we want to analyze (pre-defined in flu_data_tokens.tsv)
data_token <- "all_multiome_paired_minus_0_sample"
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
sample_metadata <- sample_metadata[sample_metadata$aliquot %in% sample_id_list,]
high_viral_load_samples <- sort(sample_metadata[sample_metadata$viral_load == "high",]$aliquot)
low_viral_load_samples <- sort(sample_metadata[sample_metadata$viral_load == "low",]$aliquot)
all_viral_load_samples <- c(high_viral_load_samples, low_viral_load_samples)
d28_samples <- sort(sample_metadata[sample_metadata$time_point == "2_D28",]$aliquot)
d_minus_1_samples <- sort(sample_metadata[sample_metadata$time_point == "2_D_minus_1",]$aliquot)
all_day_samples <- c(d28_samples, d_minus_1_samples)
male_samples <- sort(sample_metadata[sample_metadata$sex == "M",]$aliquot)
female_samples <- sort(sample_metadata[sample_metadata$sex == "F",]$aliquot)
all_sex_samples <- c(male_samples, female_samples)

# Organize sample list by high viral load followed by low viral load
sample_id_list <- sample_id_list[order(match(sample_id_list, all_viral_load_samples))]
# Sort other types of samples according to viral load order
#all_day_samples <- all_day_samples[order(match(all_day_samples, all_viral_load_samples))]
#all_sex_samples <- all_sex_samples[order(match(all_sex_samples, all_viral_load_samples))]

# Parameters for processing both RNA-seq and ATAC data
# save_progress: If you want to save your progress
save_progress <- FALSE
# Parameters for processing RNA-seq
# record_doublets: If you want to run scDblFinder and record which cells are doublets
# run_qc: If you want to plot QC metrics (and not run the rest of the pipeline)
# run_soup: If you want to remove ambient RNA using SoupX
record_doublets <- FALSE
run_qc <- FALSE
run_soup <- FALSE
# Parameters for processing ATAC-seq data
# use_rna_labels: Uses labels from the RNA data on the ATAC-data (we assume these labels are more accurate than what the GeneIntegrationMatrix finds)
use_rna_labels <- TRUE
# subset_to_rna: Subset to only those cells we kept from RNA
subset_to_rna <- TRUE
# TODO
################## ANALYSIS ##################
if(analysis_type == "RNA_seq") {
  # Step 1 - grab matrices
  if(run_soup) {
    all_sc_exp_matrices <- process_matrices_through_soup(data_path, sample_id_list)
  } else {
    all_sc_exp_matrices <- Read_h5(data_path, sample_id_list)
  }
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
    # tagged_samples <- c("3247c65ecdbfe34a", "d360f89cf9585dfe", "48ebe8475317ba95", "3c4540710e55f7b1", "fba8595c48236db8")
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
    raw_cluster_info <- capture_cluster_info(sc_obj)
    run_differential_expression_cluster(sc_obj, marker_dir)
    # Remove messy clusters
    messy_clusters <- c(6,7,13,14,15,19,28)
    idxPass <- which(Idents(sc_obj) %in% messy_clusters)
    cellsPass <- names(sc_obj$orig.ident[-idxPass])
    sc_obj.minus.messy.clusters <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
    # Override labels manually where necessary
    #sc_obj.minus.messy.clusters <- override_cluster_label(sc_obj.minus.messy.clusters, c(15), "CD16 Mono")
    #sc_obj.minus.messy.clusters <- override_cluster_label(sc_obj.minus.messy.clusters, c(34), "pDC")
    # Remove specific samples that we're not interested in
    # sc_obj.minus.messy.clusters.removed.samples <- remove_specific_samples_from_sc_obj(sc_obj.minus.messy.clusters, tagged_samples)
    # print_UMAP(sc_obj.minus.messy.clusters.removed.samples, sample_count, "tagged", plot_dir, paste0("post.clusters_all_untagged_", date, ".png"))
    # Print info about sample representation and breakdown of categories per cell type
    print(table(sc_obj.minus.messy.clusters$sample))
    print_celltype_counts(sc_obj.minus.messy.clusters)
    # Print plots for various metadata categories for final UMAP
    # Remove noisy cells
    sc_obj.minus.messy.clusters <- remove_cells_based_on_umap(sc_obj.minus.messy.clusters, -2, 0, 1, 2.5)
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "predicted_celltype_majority_vote", plot_dir, paste0("post.clusters_by_cell_type_majority_vote_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "seurat_clusters", plot_dir, paste0("post.clusters_by_cluster_num_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "predicted.id", plot_dir, paste0("post.clusters_by_cell_type_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "viral_load", plot_dir, paste0("post.clusters_by_viral_load_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "sample", plot_dir, paste0("post.clusters_by_sample_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "day", plot_dir, paste0("post.clusters_by_day_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "sex", plot_dir, paste0("post.clusters_by_sex_", date, ".png"))
    # write cells and associated cell type majority predictions to file (for ATAC-seq labeling)
    cells_for_ATAC <- data.frame("cells" = sc_obj.minus.messy.clusters$cell_name, voted_type = sc_obj.minus.messy.clusters$predicted_celltype_majority_vote)
    write.csv(cells_for_ATAC, file = paste0(analysis_dir, "rna_seq_labeled_cells_", date, ".csv"), quote = FALSE, row.names = FALSE)
    # Combine cell types for MAGICAL and other analyses that require ATAC-seq (granularity isn't as good for ATAC-seq)
    sc_obj.minus.messy.clusters <- combine_cell_types_magical(sc_obj.minus.messy.clusters)
    # Run differential expression for each cell type within each group of interest
    run_differential_expression_group(sc_obj.minus.messy.clusters, analysis_dir, "viral_load")
    run_differential_expression_group(sc_obj.minus.messy.clusters, analysis_dir, "day")
    run_differential_expression_group(sc_obj.minus.messy.clusters, analysis_dir, "sex")
    create_magical_cell_type_proportion_file(sc_obj.minus.messy.clusters, "viral_load", high_viral_load_samples, d28_samples, male_samples)
    create_magical_cell_type_proportion_file(sc_obj.minus.messy.clusters, "day", high_viral_load_samples, d28_samples, male_samples)
    create_magical_cell_type_proportion_file(sc_obj.minus.messy.clusters, "sex", high_viral_load_samples, d28_samples, male_samples)
    create_magical_cell_type_pseudobulk_file(sc_obj.minus.messy.clusters)
  }
} else if(analysis_type == "ATAC_seq") {
  # Label input files
  inputFiles <- paste0(data_path, sample_id_list, "/outs/atac_fragments.tsv.gz")
  names(inputFiles) <- all_viral_load_samples
  # Add relevant genome for ArchR
  addArchRGenome("hg38")
  # Create ArchR project from input files
  proj <- load_archR_from_input_files(inputFiles)
  # Re-load ArchR project
  proj <- loadArchRProject(path = paste0(analysis_dir, "/ArchR/"))
  # Add sample metadata
  proj <- add_sample_metadata_atac(proj, high_viral_load_samples, low_viral_load_samples,
                                   d28_samples, d_minus_1_samples, male_samples, female_samples)
  viral_load_metadata <- parse_metadata_for_samples("viral_load")
  day_metadata <- parse_metadata_for_samples("day")
  sex_metadata <- parse_metadata_for_samples("sex")
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
  proj <- plot_atac_after_filtering(proj, date)
  # Map from scRNA reference to ATAC data
  scRNA_reference <- load_rna_reference_for_atac("~/reference/")
  proj <- map_reference_to_atac(proj)
  rm(scRNA_reference)
  plot_atac_after_integration(proj, date)
  saveArchRProject(ArchRProj = proj, load = FALSE)
  # Load ArchR project 
  proj <- loadArchRProject(path = paste0(analysis_dir, "/ArchR/"))
  proj <- add_rna_labels_for_atac_data(proj, source_rna_file = "rna_seq_labeled_cells_2023-03-22.csv", use_rna_labels, subset_to_rna)
  proj <- combine_cell_types_atac(proj)
  # If we subset to RNA, we don't need to do any majority voting in clusters, etc.
  # Otherwise, we do!
  if(subset_to_rna) {
    final_proj <- proj
    # Just for convenience in code below
    final_proj$Cell_type_voting <- final_proj$predictedGroup
    plot_atac_after_majority_vote_or_subset(final_proj, date)
  } else {
    proj <- perform_majority_vote(proj)
    plot_atac_after_majority_vote_or_subset(proj, date)
    proj <- remove_single_cell_clusters(proj)
    cluster_info <- get_cluster_info(proj)
    proj <- override_cluster_label(proj, c("C8"), "CD16 Mono")
    
  }
  
  
  # If we didn't subset to RNA, then we need to do some kind of majority vote in our clusters
  if(!subset_to_rna) {
    #Remove the messy clusters (determined through visual inspection and seeing distribution of cells in each cluster)
    idxPass <- which(proj$Clusters %in% c("C1", "C29", "C30", "C31", "C38"))
    cellsPass <- proj$cellNames[-idxPass]
    proj_minus_clusters <- proj[cellsPass, ]
    saveArchRProject(ArchRProj = proj_minus_clusters, outputDirectory = paste0(analysis_dir, "ArchR_minus_clusters_rna"), load = FALSE)
    proj_minus_clusters <- loadArchRProject(path = paste0(analysis_dir, "/ArchR_minus_clusters_rna/"))
    final_proj <- proj_minus_clusters
    p1 <- plotEmbedding(ArchRProj = proj_minus_clusters, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE)
    p2 <- plotEmbedding(ArchRProj = proj_minus_clusters, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
    p3 <- plotEmbedding(ArchRProj = proj_minus_clusters, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP", force = TRUE)
    p4 <- plotEmbedding(ArchRProj = proj_minus_clusters, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", force = TRUE)
    p5 <- plotEmbedding(ArchRProj = proj_minus_clusters, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE)
    plotPDF(p1,p2,p3,p4,p5, name = paste0("Integrated_Clustering_Gene_Integration_Voting_1_RNA_labels_subset_minus_clusters_", date), ArchRProj = proj_minus_clusters, addDOC = FALSE, width = 5, height = 5)
  } else {
    final_proj <- proj
    # Just for convenience in code below
    final_proj$Cell_type_voting <- final_proj$predictedGroup
  }
  # See distribution for different metadata categories
  for (cell_type in unique(final_proj$Cell_type_voting)) {
    print(cell_type)
    idxPass <- which(final_proj$predictedGroup %in% cell_type)
    cellsPass <- final_proj$cellNames[idxPass]
    filtered_cluster <-final_proj[cellsPass,]
    print(length(cellsPass))
    print(table(filtered_cluster$viral_load))
    print(table(filtered_cluster$day))
    print(table(filtered_cluster$sex))
  }
  # Create cell type proportion file for MAGICAL
  #metadata_categories <- c("viral_load", "day", "sex")
  metadata_categories <- c("day")
  for(metadata_category in metadata_categories) {
    if(metadata_category == "day") {
      metadata <- day_metadata
    }
    cell_type_proportions_df <- data.frame("Condition" = metadata, "Sample_name" = names(metadata))
    total_cell_counts_df <- data.frame("Sample_name" = names(metadata))
    cell_counts <- vector()
    # Find total cell counts for each sample
    for (sample_id in names(metadata)) {
      idxPass <- which(final_proj$Sample %in% sample_id)
      print(length(idxPass))
      cellsPass <- final_proj$cellNames[idxPass]
      sample_subset <- subsetCells(final_proj, cellsPass)
      cell_counts <- append(cell_counts, nCells(sample_subset))
    }
    total_cell_counts_df <- cbind(total_cell_counts_df, cell_counts)
    for (cell_type in unique(final_proj$Cell_type_voting)) {
      cell_type_proportions <- vector()
      print(cell_type)
      # Grab cells associated with cell type
      idxPass <- which(final_proj$Cell_type_voting %in% cell_type)
      print(length(idxPass))
      cellsPass <- final_proj$cellNames[idxPass]
      cells_subset <- subsetCells(final_proj, cellsPass)
      for (sample_id in names(metadata)) {
        # Subset further based on cells associated with sample ID
        idxPass <- which(cells_subset$Sample %in% sample_id)
        print(length(idxPass))
        cellsPass <- cells_subset$cellNames[idxPass]
        sample_subset <- subsetCells(cells_subset, cellsPass)
        cell_counts <- nCells(sample_subset)
        cell_type_proportions <- append(cell_type_proportions, cell_counts / total_cell_counts_df[total_cell_counts_df$Sample_name == sample_id,]$cell_counts)
      }
      temp_df <- data.frame(cell_type_proportions)
      names(temp_df)[names(temp_df) == "cell_type_proportions"] <- cell_type
      cell_type_proportions_df <- cbind(cell_type_proportions_df, temp_df)
    }
    write.csv(cell_type_proportions_df, file = paste0(analysis_dir, "ATAC_cell_type_proportion_", metadata_category, ".csv"), quote = FALSE, row.names = FALSE)
    
    # Apparently info about ArchR Genome is not stored in the .RData file, so if we load an .RData file above,
    # We need to re-add our hg38 genome
    addArchRGenome("hg38")
    # Peak calling - first, call addGroupCoverages to find pseudo-bulk replicates, then call peaks using MACS2
    final_proj <- addGroupCoverages(ArchRProj = final_proj, groupBy = "Cell_type_voting", force = TRUE)
    pathToMacs2 <- findMacs2()
    final_proj <- addReproduciblePeakSet(
      ArchRProj = final_proj, 
      groupBy = "Cell_type_voting", 
      pathToMacs2 = pathToMacs2,
      force = TRUE
    )
    # Create peak matrix (matrix containing insertion counts within our merged peak set) for differential accessibility
    # calculations
    final_proj.2 <- addPeakMatrix(final_proj)
    #save.image(paste0(image_dir, "atac_after_peak_matrix.RData"))
    # Calculate differential accessible peaks for each cell type
    differential_peaks_dir <- paste0(analysis_dir, "diff_peaks/")
    if (!dir.exists(differential_peaks_dir)) {dir.create(differential_peaks_dir)}
    for (cell_type in unique(final_proj.2$Cell_type_voting)) {
      print(cell_type)
      # Grab cells associated with cell type
      idxPass <- which(final_proj.2$Cell_type_voting %in% cell_type)
      print(length(idxPass))
      cellsPass <- final_proj.2$cellNames[idxPass]
      cells_subset <- subsetCells(final_proj.2, cellsPass)
      # Find DAPs
      marker_D28_D1 <- getMarkerFeatures(ArchRProj = cells_subset, useMatrix = "PeakMatrix", groupBy = "day",
                                         testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"), maxCells = 15000,
                                         useGroups = "D28", bgdGroups = "D_MINUS_1")
      # Grab relevant stats
      marker_log_2_fc <- assays(marker_D28_D1)$Log2FC
      marker_mean <- assays(marker_D28_D1)$Mean
      marker_fdr <- assays(marker_D28_D1)$FDR
      marker_pval <- assays(marker_D28_D1)$Pval
      marker_mean_diff <- assays(marker_D28_D1)$MeanDiff
      marker_auc <- assays(marker_D28_D1)$AUC
      marker_mean_bgd <- assays(marker_D28_D1)$MeanBGD
      # Print stats to Excel spreadsheet (used by MAGICAL)
      cell_type <- sub(" ", "_", cell_type)
      list_of_datasets <- list("Log2FC" = marker_log_2_fc, "Mean" = marker_mean, "FDR" = marker_fdr, "Pval" = marker_pval, 
                               "MeanDiff" = marker_mean_diff, "AUC" = marker_auc, "MeanBGD" = marker_mean_bgd)
      write.xlsx(list_of_datasets, file = paste0(differential_peaks_dir, cell_type, "_", "HVL_LVL_diff.xlsx"), colNames = FALSE)
    }
    
    # Create Peaks.txt file for MAGICAL
    peaks <- getPeakSet(final_proj.2)
    seq_names <- as.data.frame(peaks@seqnames)
    ranges <- as.data.frame(peaks@ranges)
    peak_txt_file <- cbind(seq_names, ranges)
    colnames(peak_txt_file)[1] <- "chr"
    peak_txt_file <- peak_txt_file[, !(names(peak_txt_file) %in% c("width", "names"))]
    peak_txt_file[,1] <- sub("chr", "", peak_txt_file[,1])
    peak_txt_file[,1] <- sub("X", "23", peak_txt_file[,1])
    write.table(peak_txt_file, file = paste0(analysis_dir, "Peaks.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    
    # Create peak_motif_matches.txt file for MAGICAL
    final_proj.2 <- addMotifAnnotations(ArchRProj = final_proj.2, motifSet = "cisbp", name = "Motif")
    peak_motif_matches <- getMatches(final_proj.2, name = "Motif")
    peak_motif_txt_file <- as.data.frame(peak_motif_matches@assays@data$matches)
    # Remove _.* from ends of column names
    for (col in 1:ncol(peak_motif_txt_file)){
      colnames(peak_motif_txt_file)[col] <-  sub("_.*", "", colnames(peak_motif_txt_file)[col])
    }
    peak_motif_txt_file <- cbind(peak_txt_file, peak_motif_txt_file)
    colnames(peak_motif_txt_file)[2] <- "point1"
    colnames(peak_motif_txt_file)[3] <- "point2"
    cols <- sapply(peak_motif_txt_file, is.logical)
    peak_motif_txt_file[,cols] <- lapply(peak_motif_txt_file[,cols], as.numeric)
    write.table(peak_motif_txt_file, file = paste0(analysis_dir, "peak_motif_matches.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
    
    # Create text file for each cell type containing pseudobulk counts for peaks
    Cell_types <- unique(final_proj$Cell_type_voting)
    sample.names <- unique(final_proj$Sample)
    peak_count <- getMatrixFromProject(ArchRProj = final_proj.2, useMatrix = "PeakMatrix", useSeqnames = NULL, verbose = TRUE,binarize = FALSE,threads = getArchRThreads(),logFile = createLogFile("getMatrixFromProject"))
    pseudo_bulk_dir <- paste0(analysis_dir, "pseudo_bulk/")
    if (!dir.exists(pseudo_bulk_dir)) {dir.create(pseudo_bulk_dir)}
    for (i in c(1:length(Cell_types))){
      pseudo_bulk <- matrix(nrow = length(peaks), ncol = length(sample.names), 0)
      colnames(pseudo_bulk)<-sample.names
      rownames(pseudo_bulk)<-peaks@elementMetadata$idx
      
      for (s in c(1:length(sample.names))){
        idxMatch <- which(str_detect(peak_count$Cell_type_voting,Cell_types[i]) & str_detect(as.character(peak_count$Sample),sample.names[s]))
        if (length(idxMatch)>1){
          pseudo_bulk[,s] = Matrix::rowSums(peak_count@assays@data$PeakMatrix[,idxMatch])
        }
      }
      
      write.table(pseudo_bulk, file = paste0(pseudo_bulk_dir, "pseudo_bulk_ATAC_count_", Cell_types[i], ".txt"), quote = FALSE, sep = "\t", col.names = NA)
    }
  }
} else {
  stop("Invalid analysis type")
}

