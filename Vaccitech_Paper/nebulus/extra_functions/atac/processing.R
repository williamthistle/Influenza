# Perform dimensionality reduction, create UMAP, and find clusters in data
dimensionality_reduc <- function(proj) {
  addArchRThreads(threads = 8)
  proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2,  force = TRUE,
                          clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 30),
                          varFeatures = 25000, dimsToUse = 2:30)
  proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", force = TRUE)
  proj <- addClusters(input = proj, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 5, knnAssign = 30, maxClusters = NULL, force = TRUE)
  saveArchRProject(ArchRProj = proj, load = FALSE)
  return(proj)
}

# Create various ATAC plots after filtering is performed
plot_atac_after_filtering <- function(proj, date) {
  # UMAP plots colored by sample, cluster ID, and TSS enrichment
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
  p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
  p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
  p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "viral_load", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
  p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "day", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
  p6 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "sex", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
  plotPDF(p1,p2,p3,p4,p5,p6, name = paste0("Integrated_Clustering_snRNA_", date), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
}

# Create ATAC plot after gene integration matrix is used to predict cell types
plot_atac_after_integration <- function(proj, date) {
  pal <- paletteDiscrete(values = proj$predictedGroup)
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP", pal = pal, force = TRUE, keepAxis = TRUE)
  plotPDF(p1, name = paste0("Integrated_annotated_with_gene_integration_matrix_", date), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
}

# Create various ATAC plots after majority voting (or subsetting based on snRNA-seq cells) is performed
plot_atac_after_majority_vote_or_subset <- function(proj, date) {
  pal <- paletteDiscrete(values = proj$Cell_type_voting)
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
  p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
  p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
  p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
  plotPDF(p1,p2,p3,p4, name = paste0("Integrated_Clustering_Gene_Integration_Majority_Vote_or_Subset_", date), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
  p3 <- p3 + scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 20))
  ggsave(paste0(analysis_dir, "Integrated_Clustering_Gene_Integration_Majority_Vote_or_Subset_", date, "_more_UMAP_ticks.png"), device = "png", dpi = 300)
}

# Load our scRNA-seq reference for ATAC reference mapping
load_rna_reference_for_atac <- function(reference_dir) {
  scRNA_reference <- LoadH5Seurat(paste0(reference_dir, "multi.h5seurat"))
  # Remove certain cell types we're not interested in
  idx <- which(scRNA_reference$celltype.l2 %in% c("Doublet", "B intermediate", "CD4 CTL", "gdT", "dnT", "ILC"))
  scRNA_reference <- scRNA_reference[,-idx]
  idx <- which(scRNA_reference$celltype.l3 == "Treg Naive")
  scRNA_reference <- scRNA_reference[,-idx]
  # Rename CD4 Proliferating and CD8 Proliferating to T Proliferating
  idx <- which(scRNA_reference$celltype.l2 %in% c("CD4 Proliferating", "CD8 Proliferating"))
  scRNA_reference$celltype.l2[idx] <- "T Proliferating"
  return(scRNA_reference)
}

# Map from our reference scRNA-seq data to our ATAC data
map_reference_to_atac <- function(proj) {
  addArchRThreads(threads = 8)
  proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",#Harmony
    seRNA = scRNA_reference,
    addToArrow = FALSE,
    groupRNA = "celltype.l2",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore",
    normalization.method = "SCT",
    force = TRUE
  )
  return(proj)
}

# Add labels from our RNA-seq data to our ATAC data (only useful in true multiome situation)
# If subset_to_rna is true, then we subset to the cells found in our snRNA-seq data
# Otherwise, if we want to keep all of our snATAC-seq data, we use the cell types found in our snRNA-seq
# data where possible and otherwise use the cell types found from the reference mapping
add_rna_labels_for_atac_data <- function(proj, analysis_dir, source_rna_file, subset_to_rna) {
  curated_snRNA_seq_cells <- read.csv(paste0(analysis_dir, source_rna_file), comment.char = "")
  if(subset_to_rna) {
    idxPass <- which(proj$cellNames %in% curated_snRNA_seq_cells$cells) 
    cellsPass <- proj$cellNames[idxPass]
    proj <- proj[cellsPass, ]
    curated_snRNA_seq_cells <- curated_snRNA_seq_cells[curated_snRNA_seq_cells$cells %in% proj$cellNames,]
    curated_snRNA_seq_cells <- curated_snRNA_seq_cells[order(match(curated_snRNA_seq_cells$cells,proj$cellNames)),]
    snRNA_seq_cell_votes <- curated_snRNA_seq_cells$voted_type
    proj <- addCellColData(ArchRProj = proj, data = snRNA_seq_cell_votes, cells = proj$cellNames, name = "predictedGroup", force = TRUE)
  } else {
    predicted_snATAC_cells <- data.frame(cell_name = proj$cellNames, voted_type = proj$predictedGroup)
    for(current_row in 1:nrow(predicted_snATAC_cells)) {
      current_snATAC_cell <- predicted_snATAC_cells[current_row,]$cell_name
      if(current_snATAC_cell %in% curated_snRNA_seq_cells$cells) {
        current_voted_type <- curated_snRNA_seq_cells[curated_snRNA_seq_cells$cells == current_snATAC_cell,]$voted_type
        predicted_snATAC_cells[current_row,]$voted_type <- current_voted_type
      }
    }
    proj <- addCellColData(ArchRProj = proj, data = predicted_snATAC_cells$voted_type, cells = proj$cellNames, name = "predictedGroup", force = TRUE)
  }
  return(proj)
}

# ATAC data have limited granularity, so we combine some cell types (because we can't easily distinguish between them)
combine_cell_types_atac <- function(proj) {
  proj$oldPredictedGroup <- proj$predictedGroup
  # Combine cell types
  Cell_type_combined <- proj$predictedGroup
  idx <- grep("CD4 T", Cell_type_combined)
  Cell_type_combined[idx] <- "CD4 Memory"
  idx <- grep("CD8 T", Cell_type_combined)
  Cell_type_combined[idx] <- "CD8 Memory"
  idx <- grep("cDC", Cell_type_combined)
  Cell_type_combined[idx] <- "cDC"
  idx <- grep("Proliferating", Cell_type_combined)
  Cell_type_combined[idx] <- "Proliferating"
  idx <- grep("B", Cell_type_combined)
  Cell_type_combined[idx] <- "B"
  proj <- addCellColData(ArchRProj = proj, data = Cell_type_combined, cells = proj$cellNames, name = "predictedGroup", force = TRUE)
  #proj$predictedGroup <- replace(proj$predictedGroup, proj$predictedGroup == "CD4 Naive", "T Naive")
  #proj$predictedGroup <- replace(proj$predictedGroup, proj$predictedGroup == "CD8 Naive", "T Naive")
  proj$predictedGroup <- replace(proj$predictedGroup, proj$predictedGroup == "NK_CD56bright", "NK")
  proj$predictedGroup <- replace(proj$predictedGroup, proj$predictedGroup == "Treg", "CD4 Memory")
  return(proj)
}

# Perform majority vote within each cluster for cell type
perform_majority_vote <- function(proj) {
  # First voting scheme
  cM <- as.matrix(confusionMatrix(proj$Clusters, proj$predictedGroup))
  pre_cluster <- rownames(cM)
  max_celltype <- colnames(cM)[apply(cM, 1 , which.max)]
  Cell_type_voting <- proj$Clusters
  for (m in c(1:length(pre_cluster))){
    idxSample <- which(proj$Clusters == pre_cluster[m])
    Cell_type_voting[idxSample] <- max_celltype[m]
  }
  proj <- addCellColData(ArchRProj = proj, data = Cell_type_voting, cells = proj$cellNames, name = "Cell_type_voting", force = TRUE)
  return(proj)
}

# Remove any clusters that are only one cell (this comes up sometimes and ArchR doesn't handle them well)
remove_single_cell_clusters <- function(proj) {
  # Remove any clusters that only have 1 cell
  unique_cluster_ids <- unique(proj$Clusters)
  unique_cluster_ids <- unique_cluster_ids[order(nchar(unique_cluster_ids), unique_cluster_ids)]
  for (cluster in unique_cluster_ids) {
    idxPass <- which(proj$Clusters %in% cluster)
    if(length(idxPass) == 1) {
      cellsPass <- proj$cellNames[-idxPass]
      proj <- proj[cellsPass, ]
    }
  }
  return(proj)
}

# Get info about distributions within each cluster (cell types, viral load, day, sex)
get_cluster_info <- function(proj) {
  cluster_cell_type_predictions <- vector()
  cluster_cell_type_distributions <- list()
  cluster_sample_distributions <- list()
  cluster_viral_load_distributions <- list()
  cluster_day_distributions <- list()
  cluster_sex_distributions <- list()
  idx <- 1
  unique_cluster_ids <- unique(proj$seurat_clusters)
  unique_cluster_ids <- unique_cluster_ids[order(nchar(unique_cluster_ids), unique_cluster_ids)]
  for (cluster in unique_cluster_ids) {
    idxPass <- which(proj$seurat_clusters %in% cluster)
    cellsPass <- proj$cellNames[idxPass]
    filtered_cluster <-proj[cellsPass,]
    cluster_cell_type_predictions <- append(cluster_cell_type_predictions, table(filtered_cluster$Cell_type_voting))
    cluster_cell_type_distributions[[idx]] <- table(filtered_cluster$predictedGroup)
    cluster_sample_distributions[[idx]] <- table(filtered_cluster$Sample)
    cluster_viral_load_distributions[[idx]] <- table(filtered_cluster$viral_load)
    cluster_day_distributions[[idx]] <- table(filtered_cluster$day)
    cluster_sex_distributions[[idx]] <- table(filtered_cluster$sex)
    idx <- idx + 1
  }
  names(cluster_cell_type_predictions) <- paste(unique_cluster_ids, "-", names(cluster_cell_type_predictions))
  return(list(cluster_cell_type_predictions, cluster_cell_type_distributions, cluster_sample_distributions, cluster_viral_load_distributions, cluster_day_distributions, cluster_sex_distributions))
}

# Sometimes, we want to manually override the cell type for a given cluster (if we feel like majority vote got it wrong)
override_cluster_label_atac <- function(proj, cluster_identities, cluster_label) {
  idxPass <- which(proj$seurat_clusters %in% cluster_identities)
  proj$Cell_type_voting[idxPass] <- cluster_label
  return(proj)
}

# We may want to remove certain cell types from our data if they aren't well defined
remove_cell_types <- function(proj, cell_types) {
  idxPass <- which(proj$Cell_type_voting %in% cell_types)
  cellsPass <- proj$cellNames[-idxPass]
  proj <- proj[cellsPass, ]
  return(proj)
}

# Method to remove cells from our dataset based on the UMAP (messy cells)
remove_cells_based_on_umap_atac <- function(proj, first_x, second_x, first_y, second_y) {
  orig.umap.coords <- getEmbedding(proj)
  orig.umap.coords$cells <- rownames(orig.umap.coords)
  deleted.cells.umap.coords <- orig.umap.coords[orig.umap.coords[,1] > first_x & orig.umap.coords[,1] < second_x,]
  deleted.cells.umap.coords <- deleted.cells.umap.coords[deleted.cells.umap.coords[,2] > first_y & deleted.cells.umap.coords[,2] < second_y,]
  final.umap.coords <- orig.umap.coords[!(orig.umap.coords$cells %in% deleted.cells.umap.coords$cells),]
  cellsPass <- rownames(final.umap.coords)
  proj <- proj[cellsPass, ]
  return(proj)
}

# Method to print viral load / day / sex distributions within each cell type
print_cell_type_distributions <- function(proj) {
  # See distribution for different metadata categories
  for (cell_type in unique(proj$Cell_type_voting)) {
    print(cell_type)
    idxPass <- which(proj$predictedGroup %in% cell_type)
    cellsPass <- proj$cellNames[idxPass]
    filtered_cluster <-proj[cellsPass,]
    print(length(cellsPass))
    print(table(filtered_cluster$subject_id))
    print(table(filtered_cluster$viral_load))
    print(table(filtered_cluster$time_point))
    print(table(filtered_cluster$sex))
  }
}

# Method to create the cell type proportion file for MAGICAL
create_cell_type_proportion_MAGICAL_atac <- function(proj, analysis_dir, sample_metadata_for_SPEEDI_df) {
  cell_type_proportions_df <- data.frame("Condition" = sample_metadata_for_SPEEDI_df$time_point, 
                                         "Sample_name" = rownames(sample_metadata_for_SPEEDI_df))
  total_cell_counts_df <- data.frame("Sample_name" = rownames(sample_metadata_for_SPEEDI_df))
  cell_counts <- vector()
  # Find total cell counts for each sample
  for(sample_id in total_cell_counts_df$Sample_name) {
    idxPass <- which(proj$Sample %in% sample_id)
    print(length(idxPass))
    cellsPass <- proj$cellNames[idxPass]
    sample_subset <- subsetCells(proj, cellsPass)
    cell_counts <- append(cell_counts, nCells(sample_subset))
  }
  total_cell_counts_df <- cbind(total_cell_counts_df, cell_counts)
  for (cell_type in unique(proj$Cell_type_voting)) {
    cell_type_proportions <- vector()
    print(cell_type)
    # Grab cells associated with cell type
    idxPass <- which(proj$Cell_type_voting %in% cell_type)
    print(length(idxPass))
    cellsPass <- proj$cellNames[idxPass]
    cells_subset <- subsetCells(proj, cellsPass)
    for (sample_id in total_cell_counts_df$Sample_name) {
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
  write.csv(cell_type_proportions_df, file = paste0(analysis_dir, "ATAC_cell_type_proportion_time_point.csv"), quote = FALSE, row.names = FALSE)
}

# Method to create pseudo bulk replicates and call peaks
pseudo_bulk_replicates_and_call_peaks <- function(proj) {
  # Peak calling - first, call addGroupCoverages to find pseudo-bulk replicates, then call peaks using MACS2
  proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Cell_type_voting", force = TRUE, minCells = 50, maxCells = 500,
                            minReplicates = 3, maxReplicates = 32)
  pathToMacs2 <- findMacs2()
  proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Cell_type_voting", 
    pathToMacs2 = pathToMacs2,
    cutOff = 0.05,
    force = TRUE
  )
  return(proj)
}

# Calculate differentially accessible peaks for each cell type
calculate_daps_for_each_cell_type <- function(atac_proj, differential_peaks_dir, metadata_df) {
  pseudo_bulk_dir <- paste0(differential_peaks_dir, "associated_pseudobulk/")
  if (!dir.exists(pseudo_bulk_dir)) {dir.create(pseudo_bulk_dir, recursive = TRUE)}
  create_pseudobulk_atac(atac_proj, pseudo_bulk_dir)
  # Group used for DEGs
  group <- "time_point"
  # Calculate differential accessible peaks for each cell type
  final_lenient_de <- data.frame(Cell_Type = character(), chr = character(), start = character(), end = character(), idx = character(),
                                 sc_log2FC = character(), sc_pval = character(), sc_FDR = character(), pseudobulk_log2FC = character(),
                                 pseudobulk_pval = character()) 
  final_stricter_de <- data.frame(Cell_Type = character(), chr = character(), start = character(), end = character(), idx = character(),
                                   sc_log2FC = character(), sc_pval = character(), sc_FDR = character(), pseudobulk_log2FC = character(),
                                   pseudobulk_pval = character())  
  final_strictest_de <- data.frame(Cell_Type = character(), chr = character(), start = character(), end = character(), idx = character(),
                                 sc_log2FC = character(), sc_pval = character(), sc_FDR = character(), pseudobulk_log2FC = character(),
                                 pseudobulk_pval = character())
  current_peaks <- getPeakSet(atac_proj)
  peak_count <- getMatrixFromProject(ArchRProj = atac_proj, useMatrix = "PeakMatrix", useSeqnames = NULL, verbose = TRUE, binarize = FALSE, threads = getArchRThreads(), logFile = createLogFile("getMatrixFromProject"))
  peak_count_matrix <- peak_count@assays@data$PeakMatrix
  all_chr <- as.data.frame(current_peaks@seqnames)
  all_coord <- as.data.frame(current_peaks@ranges)
  peak_info <- cbind(all_chr, all_coord)
  for (cell_type in unique(atac_proj$Cell_type_voting)) {
    print(cell_type)
    # Grab cells associated with cell type
    idxPass <- which(atac_proj$Cell_type_voting %in% cell_type)
    print(length(idxPass))
    cellsPass <- atac_proj$cellNames[idxPass]
    cells_subset <- atac_proj[cellsPass,]
    print(cells_subset)
    # Find DAPs (single cell level)
    marker_D28_D1 <- getMarkerFeatures(ArchRProj = cells_subset, useMatrix = "PeakMatrix", groupBy = group,
                                       testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"), normBy = "ReadsInPeaks", maxCells = 15000,
                                       useGroups = "D28", bgdGroups = "D_minus_1")
    # Grab cells for current cell type in peak matrix
    idxMatch <- which(str_detect(peak_count$Cell_type_voting,cell_type))
    cell_type_specific_peak_count_matrix <- peak_count_matrix[,idxMatch]
    # Only keep rows where peaks are expressed in greater than 1% of cells in either condition
    print("Only keeping peaks expressed in 1% of cells in either condition (SC)")
    passing_peak_indices <- find_rows_with_threshold(cell_type_specific_peak_count_matrix, metadata_df, threshold = 0.01)

    # assays(marker_D28_D1)$Pval[non_passing_peak_indices,] <- 1
    marker_de <- data.frame(chr = rownames(marker_D28_D1), start = rownames(marker_D28_D1), end = rownames(marker_D28_D1), idx = rownames(marker_D28_D1), log2FC = rownames(marker_D28_D1), mean = rownames(marker_D28_D1),
                                  mean_diff = rownames(marker_D28_D1), AUC = rownames(marker_D28_D1), mean_BGD = rownames(marker_D28_D1), pval = rownames(marker_D28_D1), FDR = rownames(marker_D28_D1))  
    marker_de$chr <- as.character(rowData(marker_D28_D1)$seqnames)
    marker_de$start <- as.character(rowData(marker_D28_D1)$start)
    marker_de$end <- as.character(rowData(marker_D28_D1)$end)
    marker_de$idx <- as.character(rowData(marker_D28_D1)$idx)
    marker_de$log2FC <- assays(marker_D28_D1)$Log2FC[,1]
    marker_de$mean <- assays(marker_D28_D1)$Mean[,1]
    marker_de$mean_diff <- assays(marker_D28_D1)$MeanDiff[,1]
    marker_de$AUC <- assays(marker_D28_D1)$AUC[,1]
    marker_de$mean_BGD <- assays(marker_D28_D1)$MeanBGD[,1]
    marker_de$pval <- assays(marker_D28_D1)$Pval[,1]
    marker_de$FDR <- assays(marker_D28_D1)$FDR[,1]
    # Rearrange order of chromosomes
    custom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")
    marker_de$chr_custom_order <- factor(marker_de$chr, levels = custom_order)
    marker_de <- marker_de[order(marker_de$chr_custom_order), ]
    marker_de$chr_custom_order <- NULL
    # Set pval to 1 for any indices that aren't passing
    non_passing_peak_indices <- setdiff(as.numeric(rownames(marker_de)), passing_peak_indices)
    marker_de[non_passing_peak_indices,]$pval <- 1
    

    cell_type_for_file_name <- sub(" ", "_", cell_type)
    write.table(marker_de, paste0(differential_peaks_dir, cell_type_for_file_name, "_", "D28_D1_diff_sc_unfiltered.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
    marker_de_passing_fc <- marker_de[marker_de$log2FC < -0.1 | marker_de$log2FC > 0.1,]
    marker_de_lenient <- marker_de_passing_fc[marker_de_passing_fc$pval < 0.05,]
    write.table(marker_de_lenient, paste0(differential_peaks_dir, cell_type_for_file_name, "_", "D28_D1_diff_sc_filtered.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
    # Perform pseudobulk correction
    pseudobulk_counts <- read.table(paste0(pseudo_bulk_dir, "pseudo_bulk_ATAC_count_", cell_type_for_file_name, ".txt"), sep = "\t", header = TRUE, check.names = FALSE)
    sample_metadata <- cells_subset$Sample
    for(j in 1:nrow(metadata_df)) {
      sample_metadata <- gsub(rownames(metadata_df)[j], metadata_df[j,1], sample_metadata)
    }
    cells_subset <- addCellColData(ArchRProj = cells_subset, data = sample_metadata, cells = cells_subset$cellNames, name = "subject_id", force = TRUE)
    pseudobulk_metadata <- sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$subject_id %in% cells_subset$subject_id,]
    pseudobulk_metadata$aliquots <- rownames(pseudobulk_metadata)
    pseudobulk_metadata <- pseudobulk_metadata[match(colnames(pseudobulk_counts), pseudobulk_metadata$aliquots),]
    pseudobulk_analysis <- DESeq2::DESeqDataSetFromMatrix(countData = pseudobulk_counts, colData = pseudobulk_metadata, design = stats::formula(paste("~ subject_id + ",group)))
    pseudobulk_analysis <- DESeq2::DESeq(pseudobulk_analysis)
    pseudobulk_analysis_results_contrast <- utils::tail(DESeq2::resultsNames(pseudobulk_analysis), n=1)
    print(pseudobulk_analysis_results_contrast)
    pseudobulk_analysis_results <- DESeq2::results(pseudobulk_analysis, name=pseudobulk_analysis_results_contrast)
    pseudobulk_analysis_results$chr <- peak_info$value
    pseudobulk_analysis_results$start <- peak_info$start
    pseudobulk_analysis_results$end <- peak_info$end
    pseudobulk_analysis_results <- pseudobulk_analysis_results[,c('chr','start','end','baseMean', "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    write.table(pseudobulk_analysis_results, paste0(differential_peaks_dir, cell_type_for_file_name, "_", "D28_D1_diff_pseudo_unfiltered.tsv"), quote = FALSE, sep = "\t")
    # Only keep rows where peaks are expressed in greater than 1% of cells in either condition
    print("Only keeping peaks expressed in 1% of cells in either condition (pseudobulk)")
    pseudobulk_analysis_results$pvalue[non_passing_peak_indices] <- 1
    # Create a marker DE object that uses pseudobulk info - this will make it easier for us to make plots, run motif enrichment, etc.
    marker_D28_D1_pseudo <- marker_D28_D1
    # Rearrange pseudobulk_analysis_results so chrs are in same order as marker_D28_D1_pseudo
    custom_order <- c("chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21", "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chrX")
    pseudobulk_analysis_results$chr_custom_order <- factor(pseudobulk_analysis_results$chr, levels = custom_order)
    pseudobulk_analysis_results <- pseudobulk_analysis_results[order(pseudobulk_analysis_results$chr_custom_order), ]
    pseudobulk_analysis_results$chr_custom_order <- NULL
    # Log2FC
    pseudo_fc <- data.frame(D28 = pseudobulk_analysis_results$log2FoldChange)
    pseudo_fc$D28[is.na(pseudo_fc$D28)] <- 0
    assays(marker_D28_D1_pseudo)$Log2FC <- pseudo_fc
    # Mean
    pseudo_baseMean <- data.frame(D28 = pseudobulk_analysis_results$baseMean)
    pseudo_baseMean$D28[is.na(pseudo_baseMean$D28)] <- 0
    assays(marker_D28_D1_pseudo)$Mean <- pseudo_baseMean
    # pval
    pseudo_pval <- data.frame(D28 = pseudobulk_analysis_results$pvalue)
    pseudo_pval$D28[is.na(pseudo_pval$D28)] <- 1
    assays(marker_D28_D1_pseudo)$Pval <- pseudo_pval
    # FDR 
    pseudo_fdr <- data.frame(D28 = pseudobulk_analysis_results$padj)
    assays(marker_D28_D1_pseudo)$FDR <- pseudo_fdr
    # Filter pseudobulk
    pseudobulk_analysis_results <- pseudobulk_analysis_results[rowSums(is.na(pseudobulk_analysis_results)) == 0, ] # Remove NAs
    pseudobulk_analysis_results <- pseudobulk_analysis_results[pseudobulk_analysis_results$pvalue < 0.05,]
    write.table(pseudobulk_analysis_results, paste0(differential_peaks_dir, cell_type_for_file_name, "_", "D28_D1_diff_pseudo_filtered.tsv"), quote = FALSE, sep = "\t")
    for(current_pseudobulk_peak_row_index in 1:nrow(pseudobulk_analysis_results)) {
      current_pseudobulk_row <- pseudobulk_analysis_results[current_pseudobulk_peak_row_index,]
      chr_peak_index_combo <- rownames(current_pseudobulk_row)
      chr <- strsplit(chr_peak_index_combo, "_")[[1]][1]
      peak_index <- strsplit(chr_peak_index_combo, "_")[[1]][2]
      # PSEUDOBULK TEST
      current_row <- data.frame(cell_type, chr, peak_index)
      names(current_row) <- c("Cell_Type", "chr", "idx")
      # Check for overlap between pseudo and lenient
      marker_de_subset_lenient <- marker_de_lenient[marker_de_lenient$chr == chr,]
      marker_de_subset_lenient <- marker_de_subset_lenient[marker_de_subset_lenient$idx == peak_index,]
      if(nrow(marker_de_subset_lenient) > 0) {
        start <- marker_de_subset_lenient$start
        end <- marker_de_subset_lenient$end
        sc_log2FC <- marker_de_subset_lenient$log2FC
        sc_pval <- marker_de_subset_lenient$pval
        sc_FDR <- marker_de_subset_lenient$FDR
        pseudobulk_log2FC <- current_pseudobulk_row$log2FoldChange
        pseudobulk_pval <- current_pseudobulk_row$pvalue
        current_row <- data.frame(cell_type, chr, start, end, peak_index, sc_log2FC, sc_pval, sc_FDR, pseudobulk_log2FC, pseudobulk_pval)
        names(current_row) <- c("Cell_Type", "chr", "start", "end", "idx", "sc_log2FC", "sc_pval", "sc_FDR", "pseudobulk_log2FC", "pseudobulk_pval")
        final_lenient_de <- rbind(final_lenient_de, current_row)
        # stricter uses p-value < 0.01 (more stringent than p-value < 0.05)
        if(sc_pval < 0.01) {
          final_stricter_de <- rbind(final_stricter_de, current_row)
        }
        # strictest uses FDR < 0.1 (more stringent than unadjusted p-value)
        if(sc_FDR < 0.1) {
          final_strictest_de <- rbind(final_strictest_de, current_row)
        }
      }
    }
    cell_type_subset_de <- final_lenient_de[final_lenient_de$Cell_Type == cell_type,]
    pos_cell_type_subset_de <- cell_type_subset_de[cell_type_subset_de$sc_log2FC > 0,]
    neg_cell_type_subset_de <- cell_type_subset_de[cell_type_subset_de$sc_log2FC < 0,]
    # POSITIVE
    # STEP 1: Use all positive peaks (no pseudobulk adjustment) for motif enrichment
    # Because we don't have too many peaks, including more noise may be worth it
    motifsUp_all_positive <- peakAnnoEnrichment(
      seMarker = marker_D28_D1,
      ArchRProj = atac_proj,
      peakAnnotation = "Motif",
      cutOff = "Pval < 0.05 & Log2FC > 0.3"
    )
    df_up_all_positive <- data.frame(TF = rownames(motifsUp_all_positive), enrichment = assays(motifsUp_all_positive)$Enrichment$D28, 
                                     mlog10Padj = assays(motifsUp_all_positive)$mlog10Padj$D28)
    df_up_all_positive <- df_up_all_positive[order(df_up_all_positive$mlog10Padj, decreasing = TRUE),]
    df_up_all_positive$p_adj <- 10 ** -df_up_all_positive$mlog10Padj
    df_up_all_positive$rank <- seq_len(nrow(df_up_all_positive))
    write.table(df_up_all_positive, paste0(differential_peaks_dir, cell_type_for_file_name, "_", "D28_D1_motif_up_all_positive.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
    ggUp <- ggplot(df_up_all_positive, aes(rank, mlog10Padj, color = mlog10Padj)) + 
      geom_point(size = 1) +
      ggrepel::geom_label_repel(
        data = df_up_all_positive[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
      ) + theme_ArchR() + 
      ylab("-log10(P-adj) Motif Enrichment") + 
      xlab("Rank Sorted TFs Enriched") +
      scale_color_gradientn(colors = paletteContinuous(set = "comet"))
    ggplot2::ggsave(filename = paste0(differential_peaks_dir, cell_type_for_file_name, "_D28_D1_motif_up_all_positive.png"), 
                    plot = ggUp, device = "png", width = 4, height = 4, 
                    units = "in")
    # STEP 2: Use positive peaks with pseudobulk correction
    # Maybe the best middle ground? Get to use subject ID in model
    motifsUp <- peakAnnoEnrichment(
      seMarker = marker_D28_D1_pseudo,
      ArchRProj = atac_proj,
      peakAnnotation = "Motif",
      cutOff = "Pval < 0.05 & Log2FC > 0"
    )
    df_up <- data.frame(TF = rownames(motifsUp), enrichment = assays(motifsUp)$Enrichment$D28, 
                        mlog10Padj = assays(motifsUp)$mlog10Padj$D28)
    df_up <- df_up[order(df_up$mlog10Padj, decreasing = TRUE),]
    df_up$p_adj <- 10 ** -df_up$mlog10Padj
    df_up$rank <- seq_len(nrow(df_up))
    write.table(df_up, paste0(differential_peaks_dir, cell_type_for_file_name, "_", "D28_D1_motif_up_pseudobulk_only.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
    ggUp <- ggplot(df_up, aes(rank, mlog10Padj, color = mlog10Padj)) + 
      geom_point(size = 1) +
      ggrepel::geom_label_repel(
        data = df_up[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
      ) + theme_ArchR() + 
      ylab("-log10(P-adj) Motif Enrichment") + 
      xlab("Rank Sorted TFs Enriched") +
      scale_color_gradientn(colors = paletteContinuous(set = "comet"))
    ggplot2::ggsave(filename = paste0(differential_peaks_dir, cell_type_for_file_name, "_D28_D1_motif_up_pseudobulk_only.png"), 
                    plot = ggUp, device = "png", width = 4, height = 4, 
                    units = "in")
    # STEP 3: Use positive peaks with overlapping SC and pseudobulk
    # Better approach if we have more samples, but maybe not the best fit here
    pos_peak_indices <- c()
    # Find indices of peaks that overlap with cell type peaks
    for(current_row_idx in 1:nrow(pos_cell_type_subset_de)) {
      current_row <- pos_cell_type_subset_de[current_row_idx,]
      current_idx <- which(marker_de$chr == current_row$chr & marker_de$idx == current_row$idx)
      pos_peak_indices <- c(pos_peak_indices, current_idx)
    }
    pos_peak_indices <- sort(pos_peak_indices)
    print(length(pos_peak_indices))
    motifsUp <- peakAnnoEnrichment_mine(
      seMarker = marker_D28_D1,
      ArchRProj = atac_proj,
      peakAnnotation = "Motif",
      background = "bgdPeaks",
      idx = pos_peak_indices
    )
    df_up <- data.frame(TF = rownames(motifsUp), enrichment = assays(motifsUp)$Enrichment$D28, 
                        mlog10Padj = assays(motifsUp)$mlog10Padj$D28)
    df_up <- df_up[order(df_up$mlog10Padj, decreasing = TRUE),]
    df_up$p_adj <- 10 ** -df_up$mlog10Padj
    df_up$rank <- seq_len(nrow(df_up))
    write.table(df_up, paste0(differential_peaks_dir, cell_type_for_file_name, "_", "D28_D1_motif_up_sc_pseudobulk_overlap.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
    ggUp <- ggplot(df_up, aes(rank, mlog10Padj, color = mlog10Padj)) + 
      geom_point(size = 1) +
      ggrepel::geom_label_repel(
        data = df_up[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
      ) + theme_ArchR() + 
      ylab("-log10(P-adj) Motif Enrichment") + 
      xlab("Rank Sorted TFs Enriched") +
      scale_color_gradientn(colors = paletteContinuous(set = "comet"))
    ggplot2::ggsave(filename = paste0(differential_peaks_dir, cell_type_for_file_name, "_D28_D1_motif_up_sc_pseudobulk_overlap.png"), 
                    plot = ggUp, device = "png", width = 4, height = 4, 
                    units = "in")    
    # NEGATIVE
    # STEP 1: Use all negative peaks (no pseudobulk adjustment) for motif enrichment
    # Because we don't have too many peaks, including more noise may be worth it
    motifsDown_all_negative <- peakAnnoEnrichment(
      seMarker = marker_D28_D1,
      ArchRProj = atac_proj,
      peakAnnotation = "Motif",
      cutOff = "Pval < 0.05 & Log2FC < -0.3"
    )
    df_down_all_negative <- data.frame(TF = rownames(motifsDown_all_negative), enrichment = assays(motifsDown_all_negative)$Enrichment$D28, 
                                       mlog10Padj = assays(motifsDown_all_negative)$mlog10Padj$D28)
    df_down_all_negative <- df_down_all_negative[order(df_down_all_negative$mlog10Padj, decreasing = TRUE),]
    df_down_all_negative$p_adj <- 10 ** -df_down_all_negative$mlog10Padj
    df_down_all_negative$rank <- seq_len(nrow(df_down_all_negative))
    write.table(df_down_all_negative, paste0(differential_peaks_dir, cell_type_for_file_name, "_", "D28_D1_motif_down_all_negative.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
    ggUp <- ggplot(df_down_all_negative, aes(rank, mlog10Padj, color = mlog10Padj)) + 
      geom_point(size = 1) +
      ggrepel::geom_label_repel(
        data = df_down_all_negative[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
      ) + theme_ArchR() + 
      ylab("-log10(P-adj) Motif Enrichment") + 
      xlab("Rank Sorted TFs Enriched") +
      scale_color_gradientn(colors = paletteContinuous(set = "comet"))
    ggplot2::ggsave(filename = paste0(differential_peaks_dir, cell_type_for_file_name, "_D28_D1_motif_down_all_negative.png"), 
                    plot = ggUp, device = "png", width = 4, height = 4, 
                    units = "in")
    # STEP 2: Use negative peaks with pseudobulk correction
    # Maybe the best middle ground? Get to use subject ID in model
    motifsDown <- peakAnnoEnrichment(
      seMarker = marker_D28_D1_pseudo,
      ArchRProj = atac_proj,
      peakAnnotation = "Motif",
      background = "bgdPeaks",
      cutOff = "Pval < 0.05 & Log2FC < 0"
    )
    df_down <- data.frame(TF = rownames(motifsDown), enrichment = assays(motifsDown)$Enrichment$D28, 
                          mlog10Padj = assays(motifsDown)$mlog10Padj$D28)
    df_down <- df_down[order(df_down$mlog10Padj, decreasing = TRUE),]
    df_down$p_adj <- 10 ** -df_down$mlog10Padj
    df_down$rank <- seq_len(nrow(df_down))
    write.table(df_down, paste0(differential_peaks_dir, cell_type_for_file_name, "_", "D28_D1_motif_down_pseudobulk_only.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
    ggdown <- ggplot(df_down, aes(rank, mlog10Padj, color = mlog10Padj)) + 
      geom_point(size = 1) +
      ggrepel::geom_label_repel(
        data = df_down[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
      ) + theme_ArchR() + 
      ylab("-log10(P-adj) Motif Enrichment") + 
      xlab("Rank Sorted TFs Enriched") +
      scale_color_gradientn(colors = paletteContinuous(set = "comet"))
    ggplot2::ggsave(filename = paste0(differential_peaks_dir, cell_type_for_file_name, "_D28_D1_motif_down_pseudobulk_only.png"), 
                    plot = ggdown, device = "png", width = 4, height = 4, 
                    units = "in")
    # STEP 3: Use negative peaks with overlapping SC and pseudobulk
    # Better approach if we have more samples, but maybe not the best fit here
    neg_peak_indices <- c()
    # Find indices of peaks that overlap with cell type peaks
    for(current_row_idx in 1:nrow(neg_cell_type_subset_de)) {
      current_row <- neg_cell_type_subset_de[current_row_idx,]
      current_idx <- which(marker_de$chr == current_row$chr & marker_de$idx == current_row$idx)
      neg_peak_indices <- c(neg_peak_indices, current_idx)
    }
    neg_peak_indices <- sort(neg_peak_indices)
    print(length(neg_peak_indices))
    motifsDown <- peakAnnoEnrichment_mine(
      seMarker = marker_D28_D1,
      ArchRProj = atac_proj,
      peakAnnotation = "Motif",
      background = "bgdPeaks",
      idx = neg_peak_indices
    )
    df_down <- data.frame(TF = rownames(motifsDown), enrichment = assays(motifsDown)$Enrichment$D28, 
                          mlog10Padj = assays(motifsDown)$mlog10Padj$D28)
    df_down <- df_down[order(df_down$mlog10Padj, decreasing = TRUE),]
    df_down$p_adj <- 10 ** -df_down$mlog10Padj
    df_down$rank <- seq_len(nrow(df_down))
    write.table(df_down, paste0(differential_peaks_dir, cell_type_for_file_name, "_", "D28_D1_motif_down_sc_pseudobulk_overlap.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
    ggdown <- ggplot(df_down, aes(rank, mlog10Padj, color = mlog10Padj)) + 
      geom_point(size = 1) +
      ggrepel::geom_label_repel(
        data = df_down[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
      ) + theme_ArchR() + 
      ylab("-log10(P-adj) Motif Enrichment") + 
      xlab("Rank Sorted TFs Enriched") +
      scale_color_gradientn(colors = paletteContinuous(set = "comet"))
    ggplot2::ggsave(filename = paste0(differential_peaks_dir, cell_type_for_file_name, "_D28_D1_motif_down_sc_pseudobulk_overlap.png"), 
                    plot = ggdown, device = "png", width = 4, height = 4, 
                    units = "in")    
  }
  # Write final lenient table
  write.table(final_lenient_de, paste0(differential_peaks_dir, "D28_D1_diff_lenient_final.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  pos_final_lenient_de <- final_lenient_de[final_lenient_de$sc_log2FC > 0,]
  write.table(pos_final_lenient_de, paste0(differential_peaks_dir, "D28_D1_diff_lenient_final_pos.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  neg_final_lenient_de <- final_lenient_de[final_lenient_de$sc_log2FC < 0,]
  write.table(neg_final_lenient_de, paste0(differential_peaks_dir, "D28_D1_diff_lenient_final_neg.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  # Write final stricter table
  write.table(final_stricter_de, paste0(differential_peaks_dir, "D28_D1_diff_stricter_final.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  pos_final_stricter_de <- final_stricter_de[final_stricter_de$sc_log2FC > 0,]
  write.table(pos_final_stricter_de, paste0(differential_peaks_dir, "D28_D1_diff_stricter_final_pos.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  neg_final_stricter_de <- final_stricter_de[final_stricter_de$sc_log2FC < 0,]
  write.table(neg_final_stricter_de, paste0(differential_peaks_dir, "D28_D1_diff_stricter_final_neg.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  # Write final strictest table
  write.table(final_strictest_de, paste0(differential_peaks_dir, "D28_D1_diff_strictest_final.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  pos_final_strictest_de <- final_strictest_de[final_strictest_de$sc_log2FC > 0,]
  write.table(pos_final_strictest_de, paste0(differential_peaks_dir, "D28_D1_diff_strictest_final_pos.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  neg_final_strictest_de <- final_strictest_de[final_strictest_de$sc_log2FC < 0,]
  write.table(neg_final_strictest_de, paste0(differential_peaks_dir, "D28_D1_diff_strictest_final_neg.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  print(paste0("Done performing differential expression for group ", group, " for each cell type"))
}

# Create Peaks.txt file (list of peaks) used by MAGICAL
create_peaks_file <- function(proj, analysis_dir) {
  peaks <- getPeakSet(proj)
  seq_names <- as.data.frame(peaks@seqnames)
  ranges <- as.data.frame(peaks@ranges)
  peak_txt_file <- cbind(seq_names, ranges)
  colnames(peak_txt_file)[1] <- "chr"
  peak_txt_file <- peak_txt_file[, !(names(peak_txt_file) %in% c("width", "names"))]
  peak_txt_file[,1] <- sub("chr", "", peak_txt_file[,1])
  peak_txt_file[,1] <- sub("X", "23", peak_txt_file[,1])
  write.table(peak_txt_file, file = paste0(analysis_dir, "Peaks.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
  return(peak_txt_file)
}

# Create peak_motif_matches.txt file for MAGICAL
create_peak_motif_matches_file <- function(proj, analysis_dir, peak_txt_file) {
  proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
  peak_motif_matches <- getMatches(proj, name = "Motif")
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
  return(proj)
}

# Create pseudobulk count files for each cell type
create_pseudobulk_atac <- function(proj, pseudo_bulk_dir) {
  # Create text file for each cell type containing pseudobulk counts for peaks
  Cell_types <- unique(proj$Cell_type_voting)
  sample.names <- unique(proj$Sample)
  peaks <- getPeakSet(proj)
  peak_count <- getMatrixFromProject(ArchRProj = proj, useMatrix = "PeakMatrix", useSeqnames = NULL, verbose = TRUE, binarize = FALSE, threads = getArchRThreads(), logFile = createLogFile("getMatrixFromProject"))
  for (i in c(1:length(Cell_types))){
    pseudo_bulk <- matrix(nrow = length(peaks), ncol = length(sample.names), 0)
    # We add Sample_ to the beginning of each column name to avoid MATLAB (MAGICAL) complaining - currently not doing this
    colnames(pseudo_bulk)<- sample.names
    rownames(pseudo_bulk) <- paste0(as.character(peaks@seqnames), "_", peaks@elementMetadata$idx)
    for (s in c(1:length(sample.names))){
      idxMatch <- which(str_detect(peak_count$Cell_type_voting,Cell_types[i]) & str_detect(as.character(peak_count$Sample),sample.names[s]))
      if (length(idxMatch)>1){
        pseudo_bulk[,s] = Matrix::rowSums(peak_count@assays@data$PeakMatrix[,idxMatch])
      }
    }
    cell_type <- sub(" ", "_", Cell_types[i])
    write.table(pseudo_bulk, file = paste0(pseudo_bulk_dir, "pseudo_bulk_ATAC_count_", cell_type, ".txt"), quote = FALSE, sep = "\t")
  }
}

# Create pseudobulk count files for each cell type - new? approach
create_pseudobulk_atac_new <- function(proj, metadata_df, pseudo_bulk_dir) {
  # Create text file for each cell type containing pseudobulk counts for peaks
  Cell_types <- unique(proj$Cell_type_voting)
  sample.names <- unique(proj$Sample)
  peaks <- getPeakSet(proj)
  peak_count <- getMatrixFromProject(ArchRProj = proj, useMatrix = "PeakMatrix", useSeqnames = NULL, verbose = TRUE, binarize = FALSE, threads = getArchRThreads(), logFile = createLogFile("getMatrixFromProject"))
  peak_count_matrix <- peak_count@assays@data$PeakMatrix
  for (i in c(1:length(Cell_types))){
    print(Cell_types[i])
    # Grab cells for current cell type in peak matrix
    print("Grabbing cells for current cell type")
    idxMatch <- which(str_detect(peak_count$Cell_type_voting,Cell_types[i]))
    cell_type_specific_peak_count_matrix <- peak_count_matrix[,idxMatch]
    # Only keep rows where peaks are expressed in greater than 5% of cells in either condition
    print("Only keeping peaks expressed in 5% of cells in either condition")
    passing_peak_indices <- find_rows_with_threshold(cell_type_specific_peak_count_matrix, metadata_df, threshold = 0.05)
    cell_type_specific_peak_count_matrix <- cell_type_specific_peak_count_matrix[passing_peak_indices,]
    pseudo_bulk <- matrix(nrow = nrow(cell_type_specific_peak_count_matrix), ncol = length(sample.names), 0)
    colnames(pseudo_bulk) <- sample.names
    cell_type_peaks <- peaks[passing_peak_indices,]
    rownames(pseudo_bulk) <- paste0(as.character(cell_type_peaks@seqnames), "_", cell_type_peaks@elementMetadata$idx)
    print("Adding peak counts for each sample to pseudobulk matrix")
    for (s in c(1:length(sample.names))){
      print(sample.names[s])
      idxMatch <- grep(sample.names[s], colnames(cell_type_specific_peak_count_matrix))
      if (length(idxMatch)>1){
        pseudo_bulk[,s] <- Matrix::rowSums(cell_type_specific_peak_count_matrix[,idxMatch])
      } else {
        pseudo_bulk[,s] <- 0
      }
    }
    cell_type <- sub(" ", "_", Cell_types[i])
    write.table(pseudo_bulk, file = paste0(pseudo_bulk_dir, "pseudo_bulk_ATAC_count_", cell_type, ".txt"), quote = FALSE, sep = "\t")
  }
}






# Create pseudobulk count files (sample level)
create_pseudobulk_atac_sample <- function(proj) {
  sample.names <- unique(proj$Sample)
  peaks <- getPeakSet(proj)
  peak_count <- getMatrixFromProject(ArchRProj = proj, useMatrix = "PeakMatrix", useSeqnames = NULL, verbose = TRUE, binarize = FALSE, threads = getArchRThreads(), logFile = createLogFile("getMatrixFromProject"))
  pseudo_bulk <- matrix(nrow = length(peaks), ncol = length(sample.names), 0)
  colnames(pseudo_bulk)<- sample.names
  rownames(pseudo_bulk) <- paste0(as.character(peaks@seqnames), "_", peaks@elementMetadata$idx)
  for (s in c(1:length(sample.names))){
    idxMatch <- which(str_detect(as.character(peak_count$Sample),sample.names[s]))
    if (length(idxMatch)>1){
      pseudo_bulk[,s] = Matrix::rowSums(peak_count@assays@data$PeakMatrix[,idxMatch])
    }
  }
  return(pseudo_bulk)
}


create_magical_input_files_atac <- function(atac_proj, MAGICAL_file_dir) {
  if (!dir.exists(MAGICAL_file_dir)) {dir.create(MAGICAL_file_dir)}
  MAGICAL_cell_metadata_dir <- paste0(MAGICAL_file_dir, "scATAC_Cell_Metadata/")
  if (!dir.exists(MAGICAL_cell_metadata_dir)) {dir.create(MAGICAL_cell_metadata_dir)}
  MAGICAL_peak_counts_dir <- paste0(MAGICAL_file_dir, "scATAC_Read_Counts/")
  if (!dir.exists(MAGICAL_peak_counts_dir)) {dir.create(MAGICAL_peak_counts_dir)}
  MAGICAL_candidate_peaks_dir <- paste0(MAGICAL_file_dir, "scATAC_Peak_Coordinates/")
  if (!dir.exists(MAGICAL_candidate_peaks_dir)) {dir.create(MAGICAL_candidate_peaks_dir)}
  MAGICAL_motif_mapping_prior_dir <- paste0(MAGICAL_file_dir, "scATAC_Motif_Mapping_Prior/")
  if (!dir.exists(MAGICAL_motif_mapping_prior_dir)) {dir.create(MAGICAL_motif_mapping_prior_dir)}

  for(cell_type in unique(atac_proj$Cell_type_voting)) {
    print(cell_type)
    cell_type_for_file_name <- sub(" ", "_", cell_type)
    # Grab cells associated with cell type
    idxPass <- which(atac_proj$Cell_type_voting %in% cell_type)
    cellsPass <- atac_proj$cellNames[idxPass]
    # Subset ArchR project and peak matrix to associated cells
    atac_proj_cell_type_subset <- atac_proj[cellsPass,]
    # 1) Cell metadata
    atac_proj_metadata_df <- data.frame(cell_index = seq(length(atac_proj_cell_type_subset$cellNames)), 
                                       cell_barcode = atac_proj_cell_type_subset$cellNames, 
                                       cell_type = atac_proj_cell_type_subset$Cell_type_voting, 
                                       sample = atac_proj_cell_type_subset$Sample, 
                                       condition = atac_proj_cell_type_subset$time_point)
    write.table(atac_proj_metadata_df, file = paste0(MAGICAL_cell_metadata_dir, cell_type_for_file_name, "_ATAC_cell_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    # 2) ATAC assay cell count
    atac_peak_matrix_cell_type_subset <- getMatrixFromProject(ArchRProj = atac_proj_cell_type_subset, useMatrix = "PeakMatrix", useSeqnames = NULL,
                                                             verbose = TRUE,binarize = FALSE,threads = getArchRThreads(),
                                                             logFile = createLogFile("getMatrixFromProject"))
    
    final_peak_matrix <- atac_peak_matrix_cell_type_subset@assays@data@listData$PeakMatrix
    final_peak_matrix <- final_peak_matrix[, atac_proj_cell_type_subset$cellNames]
    saveRDS(final_peak_matrix, file= paste0(MAGICAL_peak_counts_dir, cell_type_for_file_name, "_ATAC_read_counts.rds"))
  }
  
  # 2) Peak set
  current_peaks <- getPeakSet(atac_proj)
  write.table(current_peaks[,1:13], file = paste0(MAGICAL_candidate_peaks_dir, "ATAC_peak_coordinates.tsv"), quote = FALSE, col.names = FALSE,  sep = "\t")
  
  atac_proj <- addMotifAnnotations(ArchRProj = atac_proj, motifSet = "cisbp", name = "Motif",
                                                  force = TRUE)
  # A Motif-Matches-In-Peaks.rds file will be created under the Annotations folder
  peak_motif_mapping <- readRDS(file = paste0(ATAC_output_dir, "minus_clusters_no_batch_correction/Annotations/Motif-Matches-In-Peaks.rds"))
  write.table(lapply(summary(peak_motif_mapping@assays@data@listData$matches), as.numeric), 
              file=paste0(MAGICAL_motif_mapping_prior_dir, "ATAC_motif_mapping_prior.tsv"),
              quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = "\t")
}

  

peakAnnoEnrichment_mine <- function(seMarker = NULL, ArchRProj = NULL, peakAnnotation = NULL, 
                                    matches = NULL, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", background = "all", idx = NULL, 
                                    logFile = createLogFile("peakAnnoEnrichment")) {
  ArchR:::.validInput(input = seMarker, name = "seMarker", valid = c("SummarizedExperiment"))
  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchR:::.validInput(input = peakAnnotation, name = "peakAnnotation", 
                      valid = c("character", "null"))
  ArchR:::.validInput(input = matches, name = "matches", valid = c("SummarizedExperiment", 
                                                                   "null"))
  ArchR:::.validInput(input = cutOff, name = "cutOff", valid = c("character"))
  ArchR:::.validInput(input = background, name = "background", valid = c("character"))
  ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))
  tstart <- Sys.time()
  ArchR:::.startLogging(logFile = logFile)
  ArchR:::.logThis(mget(names(formals()), sys.frame(sys.nframe())), 
                   "peakAnnoEnrichment Input-Parameters", logFile = logFile)
  if (metadata(seMarker)$Params$useMatrix != "PeakMatrix") {
    stop("Only markers identified from PeakMatrix can be used!")
  }
  if (is.null(matches)) {
    matches <- ArchR::getMatches(ArchRProj, peakAnnotation)
  }
  r1 <- SummarizedExperiment::rowRanges(matches)
  pr1 <- paste(seqnames(r1), start(r1), end(r1), sep = "_")
  mcols(r1) <- NULL
  r2 <- ArchR::getPeakSet(ArchRProj)
  pr2 <- paste(seqnames(r2), start(r2), end(r2), sep = "_")
  mcols(r2) <- NULL
  r3 <- GRanges(rowData(seMarker)$seqnames, IRanges(rowData(seMarker)$start, 
                                                    rowData(seMarker)$end))
  pr3 <- paste(seqnames(r3), start(r3), end(r3), sep = "_")
  mcols(r3) <- NULL
  ArchR:::.logThis(r1, "Peaks-Matches", logFile = logFile)
  ArchR:::.logThis(r2, "Peaks-ArchRProj", logFile = logFile)
  ArchR:::.logThis(r3, "Peaks-SeMarker", logFile = logFile)
  ArchR:::.logThis(pr1, "Peaks-Pasted-Matches", logFile = logFile)
  ArchR:::.logThis(pr2, "Peaks-Pasted-ArchRProj", logFile = logFile)
  ArchR:::.logThis(pr3, "Peaks-Pasted-SeMarker", logFile = logFile)
  if (length(which(pr1 %ni% pr2)) != 0) {
    stop("Peaks from matches do not match peakSet in ArchRProj!")
  }
  if (length(which(pr2 %ni% pr3)) != 0) {
    stop("Peaks from seMarker do not match peakSet in ArchRProj!")
  }
  rownames(matches) <- pr1
  matches <- matches[pr3, ]
  assayNames <- names(SummarizedExperiment::assays(seMarker))
  for (an in assayNames) {
    eval(parse(text = paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", 
                             an, "']]")))
  }
  passMat <- eval(parse(text = cutOff))
  passMat[is.na(passMat)] <- FALSE
  for (an in assayNames) {
    eval(parse(text = paste0("rm(", an, ")")))
  }
  if (tolower(background) %in% c("backgroundpeaks", "bgdpeaks", 
                                 "background", "bgd")) {
    method <- "bgd"
    bgdPeaks <- SummarizedExperiment::assay(getBgdPeaks(ArchRProj))
  }
  else {
    method <- "all"
  }
  enrichList <- lapply(seq_len(ncol(seMarker)), function(x) {
    ArchR:::.logDiffTime(sprintf("Computing Enrichments %s of %s", 
                                 x, ncol(seMarker)), t1 = tstart, verbose = TRUE, 
                         logFile = logFile)
    if(is.null(idx)) {
      idx <- which(passMat[, x])
    } else {
      names(idx) <- idx
    }
    print(idx)
    print(length(idx))
    if (method == "bgd") {
      ArchR:::.computeEnrichment(matches, idx, c(idx, as.vector(bgdPeaks[idx, 
      ])))
    }
    else {
      ArchR:::.computeEnrichment(matches, idx, seq_len(nrow(matches)))
    }
  }) %>% SimpleList
  names(enrichList) <- colnames(seMarker)
  assays <- lapply(seq_len(ncol(enrichList[[1]])), function(x) {
    d <- lapply(seq_along(enrichList), function(y) {
      enrichList[[y]][colnames(matches), x, drop = FALSE]
    }) %>% Reduce("cbind", .)
    colnames(d) <- names(enrichList)
    d
  }) %>% SimpleList
  names(assays) <- colnames(enrichList[[1]])
  assays <- rev(assays)
  out <- SummarizedExperiment::SummarizedExperiment(assays = assays)
  ArchR:::.endLogging(logFile = logFile)
  out
}

find_rows_with_threshold <- function(matrix_data, metadata_df, threshold = 0.05) {
  d28_aliquots <- rownames(metadata_df[metadata_df$time_point == "D28",])
  d28_pattern <- paste(d28_aliquots, collapse = "|")
  d_minus_1_aliquots <-  rownames(metadata_df[metadata_df$time_point == "D_minus_1",])
  d_minus_1_pattern <- paste(d_minus_1_aliquots, collapse = "|")
  d28_matrix <- matrix_data[,grep(d28_pattern, colnames(matrix_data))]
  d_minus_1_matrix <- matrix_data[grep(d_minus_1_pattern, colnames(matrix_data)),]
  row_indices_1 <- which(rowSums(d28_matrix != 0) / ncol(d28_matrix) >= threshold)
  row_indices_2 <- which(rowSums(d_minus_1_matrix != 0) / ncol(d_minus_1_matrix) >= threshold)
  row_indices <- sort(unique(row_indices_1, row_indices_2))
  return(row_indices)
}

# Run differential expression for Seurat ATAC object 
run_differential_expression_controlling_for_subject_id_atac <- function(sc_obj, analysis_dir, sample_metadata_for_SPEEDI_df, group) {
  print(paste0("Performing differential expression for group ", group, " for each cell type (controlling for subject ID)"))
  
  if (!dir.exists(analysis_dir)) {dir.create(analysis_dir, recursive = TRUE)}
  
  expected_num_samples <- length(unique(sc_obj$Sample))
  
  cell_types <- unique(sc_obj$Cell_type_voting)
  
  if(group == "viral_load") {
    first_group <- "high"
    second_group <- "low"
  } else if(group == "time_point") {
    first_group <- "D28"
    second_group <- "D_minus_1"
  } else if(group == "sex") {
    first_group <- "M"
    second_group <- "F"
  }
  
  # Make reading data parallel
  if (Sys.getenv("SLURM_NTASKS_PER_NODE") == "") {
    n.cores <- as.numeric(detectCores())
  } else {
    n.cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
  }
  
  n.cores <- 16
  
  message(paste0("Number of cores: ", n.cores))
  
  registerDoMC(n.cores)
  message("Begin parallelizing...")
  
  compiled_output <- foreach(
    i = 1:length(cell_types),
    .packages = c("Seurat", "base")
  ) %dopar% {
    current_cell_type <- cell_types[[i]]
    print(current_cell_type)
    cell_type_for_file_name <- sub(" ", "_", current_cell_type)
    # Grab indices for current cell type
    idxPass <- which(sc_obj$predicted_celltype_majority_vote %in% current_cell_type)
    if(length(idxPass) == 0) {
      print("No cells found, so skipping cell type")
    } else {
      # Grab Seurat object that only contains cells from current cell type
      cellsPass <- names(sc_obj$orig.ident[idxPass])
      cells_subset <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
      DefaultAssay(cells_subset) <- "peaks"
      # Set up differential expression
      Idents(cells_subset) <- group
      if(group == "viral_load") {
        first_group <- "high"
        second_group <- "low"
      } else if(group == "time_point") {
        first_group <- "D28"
        second_group <- "D_minus_1"
      } else if(group == "sex") {
        first_group <- "M"
        second_group <- "F"
      }
      # Perform SC DE
      current_de <- FindMarkers(cells_subset, test.use="LR", latent.vars = c('nCount_peaks', 'subject_id'), ident.1 = first_group, ident.2 = second_group, logfc.threshold = 0, min.pct = 0)
      # Save unfiltered DE results
      write.table(current_de, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type_for_file_name, "-", group, "-controlling_for_subject_id_sc_unfiltered.tsv"), quote = FALSE, sep = "\t")
      # Filter DE results and write to file
      current_de_filtered <- current_de[current_de$p_val < 0.05,]
      current_de_filtered <- current_de_filtered[abs(current_de_filtered$avg_log2FC) >= 0.1,]
      current_de_filtered <- current_de_filtered[current_de_filtered$pct.1 >= 0.01 | current_de_filtered$pct.2 >= 0.01, ]
      write.table(current_de_filtered, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type_for_file_name, "-", group, "-controlling_for_subject_id_sc_pct_0.01.tsv"), quote = FALSE, sep = "\t")
      # Run pseudobulk DE
      # Read in pseudobulk counts
      pseudobulk_counts <- create_pseudobulk_counts_atac_seurat(cells_subset)
      # Grab associated metadata
      pseudobulk_metadata <- sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$subject_id %in% cells_subset$subject_id,]
      pseudobulk_metadata$aliquots <- rownames(pseudobulk_metadata)
      pseudobulk_metadata <- pseudobulk_metadata[match(colnames(pseudobulk_counts), pseudobulk_metadata$aliquots),]
      # Set up flag matrix
      pseudobulk_counts_flag <- pseudobulk_counts 
      pseudobulk_counts_flag[pseudobulk_counts_flag == 1] <- 0
      pseudobulk_counts_flag[pseudobulk_counts_flag > 0] <- 1
      pseudobulk_counts_flag <- pseudobulk_counts_flag[rowSums(pseudobulk_counts_flag != 0) >= 2, ]
      # Scale counts
      total_ATAC_reads=5e6
      ATAC_scale_factor=total_ATAC_reads/(colSums(pseudobulk_counts)+1)
      ATAC_count=sweep(pseudobulk_counts, 2, ATAC_scale_factor, "*")
      # Log scale
      ATAC_log2=log2(ATAC_count+1)
      # Subset to relevant peaks for pseudobulk testing in pseudobulk counts
      ATAC_log2_candidate <- ATAC_log2[rownames(ATAC_log2) %in% rownames(current_de_filtered),]
      ATAC_log2_candidate <- ATAC_log2_candidate[rownames(ATAC_log2_candidate) %in% rownames(pseudobulk_counts_flag),]
      # Perform pseudobulk testing
      pseudobulk_p_values <- c()
      pseudobulk_robust_p_values <- c()
      pseudobulk_fc <- c() 
      for(current_peak_idx in 1:nrow(ATAC_log2_candidate)) {
        pseudobulk_fc <- c(pseudobulk_fc, mean(as.numeric(ATAC_log2_candidate[current_peak_idx, pseudobulk_metadata$time_point == "D28"])) - mean(as.numeric(ATAC_log2_candidate[current_peak_idx, pseudobulk_metadata$time_point == "D_minus_1"])))
        # Normal LM
        mdl <- lm(as.numeric(ATAC_log2_candidate[current_peak_idx, ]) ~ pseudobulk_metadata$time_point)
        pseudobulk_p_values <- c(pseudobulk_p_values, summary(mdl)$coefficients[2, 4])
        # Robust LM
        rsl <- MASS::rlm(as.numeric(ATAC_log2_candidate[current_peak_idx, ]) ~ pseudobulk_metadata$time_point)
        robust_pvalue <- sfsmisc::f.robftest(rsl, var = "pseudobulk_metadata$time_pointD28")$p.value
        pseudobulk_robust_p_values <- c(pseudobulk_robust_p_values, robust_pvalue)
      }
      # Create unfiltered DE results
      current_de_pseudobulk <- data.frame(logFC_pseudobulk = pseudobulk_fc, p_value_pseudobulk = pseudobulk_p_values,
                                          robust_p_value_pseudobulk = pseudobulk_robust_p_values)
      rownames(current_de_pseudobulk) <- rownames(ATAC_log2_candidate)
      # Print unfiltered DE results
      write.table(current_de_pseudobulk, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type_for_file_name, "-", group, "-controlling_for_subject_id_pseudobulk_unfiltered.tsv"), quote = FALSE, sep = "\t")
      # Filter DE results (p-value < 0.05) and write to file
      current_de_pseudobulk <- current_de_pseudobulk[rowSums(is.na(current_de_pseudobulk)) == 0, ] # Remove NAs
      current_de_pseudobulk <- current_de_pseudobulk[current_de_pseudobulk$p_value_pseudobulk < 0.05 | current_de_pseudobulk$robust_p_value_pseudobulk < 0.05,]
      write.table(current_de_pseudobulk, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type_for_file_name, "-", group, "-controlling_for_subject_id_pseudobulk.tsv"), quote = FALSE, sep = "\t")
      # We want to look at pct = 0.05 and 0.1 as well, just to see whether signal is better
      pct_thresholds <- c(0.01, 0.05, 0.1)
      for(pct_threshold in pct_thresholds) {
        # This DF will store those peaks which were found to be significant in both SC and pseudobulk
        overlapping_peak_de <- data.frame(Cell_Type = character(), Peak_Name = character(), sc_pval = character(), sc_log2FC = character(), pseudo_bulk_pval = character(), pseudo_bulk_robust_pval = character(),
                                          pseudo_bulk_log2FC = character())
        if(pct_threshold != 0.01) {
          current_de_with_current_pct <- current_de_filtered[current_de_filtered$pct.1 >= pct_threshold | current_de_filtered$pct.2 >= pct_threshold,]
          write.table(current_de_with_current_pct, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type_for_file_name, "-", group, "-controlling_for_subject_id_sc_pct_", pct_threshold, ".tsv"), quote = FALSE, sep = "\t")
        } else {
          current_de_with_current_pct <- current_de_filtered
        }
        final_peaks <- intersect(rownames(current_de_with_current_pct), rownames(current_de_pseudobulk))
        # Record information about remaining genes in overlapping_peak_de
        for(current_peak in final_peaks) {
          current_sc_pval <- current_de[rownames(current_de) == current_peak,]$p_val
          current_sc_log2FC <- current_de[rownames(current_de) == current_peak,]$avg_log2FC
          current_pseudo_bulk_pval <- current_de_pseudobulk[rownames(current_de_pseudobulk) == current_peak,]$p_value_pseudobulk
          current_pseudo_bulk_robust_pval <- current_de_pseudobulk[rownames(current_de_pseudobulk) == current_peak,]$robust_p_value_pseudobulk
          current_pseudo_bulk_log2FC <- current_de_pseudobulk[rownames(current_de_pseudobulk) == current_peak,]$logFC_pseudobulk
          current_row <- data.frame(current_cell_type, current_peak, current_sc_pval, current_sc_log2FC, current_pseudo_bulk_pval, current_pseudo_bulk_robust_pval, current_pseudo_bulk_log2FC)
          names(current_row) <- c("Cell_Type", "Peak_Name", "sc_pval", "sc_log2FC", "pseudo_bulk_pval", "pseudo_bulk_robust_pval", "pseudo_bulk_log2FC")
          overlapping_peak_de <- rbind(overlapping_peak_de, current_row)
        }
        # Make sure that FC sign matches between SC and pseudo analysis
        overlapping_peak_de <- overlapping_peak_de[sign(overlapping_peak_de$sc_log2FC) == sign(overlapping_peak_de$pseudo_bulk_log2FC),]
        # Write both positive and negative peaks to file
        write.table(overlapping_peak_de, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type_for_file_name, "-", group, "-controlling_for_subject_id_final_pct_", pct_threshold, ".tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
      }
    }
  }
  print(paste0("Done performing differential expression for group ", group, " for each cell type"))
}

create_pseudobulk_counts_atac_seurat <- function(sc_obj) {
  cells_pseudobulk <- list()
  for (sample_name in unique(sc_obj$Sample)) {
    idxMatch <- which(stringr::str_detect(as.character(sc_obj$Sample), sample_name))
    # Note - ideally, this should be >= 1, but there's a bug with Seurat V5 where data from objects with 1 cell cannot be
    # sampled correctly. Thus, in this edge case, we assume object has 0 cells
    if(length(idxMatch)>1) {
      samples_subset <- subset(x = sc_obj, subset = Sample %in% sample_name)
      samples_data <- samples_subset[["peaks"]]$data
      samples_data <- rowSums(as.matrix(samples_data))
      cells_pseudobulk[[sample_name]] <- samples_data
    } else {
      cells_pseudobulk[[sample_name]] <- numeric(nrow(sc_obj@assays$peaks))
    }
  }
  final_cells_pseudobulk_df <- dplyr::bind_cols(cells_pseudobulk[1])
  for (idx in 2:length(unique(sc_obj$Sample))) {
    final_cells_pseudobulk_df <- dplyr::bind_cols(final_cells_pseudobulk_df, cells_pseudobulk[idx])
  }
  final_cells_pseudobulk_df <- as.data.frame(final_cells_pseudobulk_df)
  rownames(final_cells_pseudobulk_df) <- names(cells_pseudobulk[[1]])
  return(final_cells_pseudobulk_df)
}

# Generate motif tables using Signac
generate_motifs_with_signac <- function(seurat_atac, motif_input_dir, motif_output_dir) {
  set.seed(SPEEDI::get_speedi_seed())
  # Declare cell types
  cell_types <- unique(seurat_atac$predicted_celltype_majority_vote)
  Idents(seurat_atac) <- "predicted_celltype_majority_vote"
  for(cell_type in cell_types) {
    print(cell_type)
    # Create dir for cell type results
    cell_type_for_file_name <- sub(" ", "_", cell_type)
    cell_type_dir <- paste0(motif_output_dir, cell_type_for_file_name, "/")
    if (!dir.exists(cell_type_dir)) {dir.create(cell_type_dir, recursive = TRUE)}
    # Grab cells associated with cell type
    idxPass <- which(seurat_atac$predicted_celltype_majority_vote %in% cell_type)
    cellsPass <- names(seurat_atac$orig.ident[idxPass])
    cells_subset <- subset(x = seurat_atac, subset = cell_name %in% cellsPass)
    DefaultAssay(cells_subset) <- "peaks"
    
    # We will use peaks open in cell type of interest for background peaks
    open.peaks <- AccessiblePeaks(seurat_atac, idents = cell_type)
    meta.feature <- GetAssayData(cells_subset, assay = "peaks", slot = "meta.features")
    
    # Analysis types include those peaks that passed SC filtering stage as well as those peaks that passed both SC filtering and pseudobulk
    analysis_types <- c("sc", "final")
    # Tested pct levels include 0.01, 0.05, 0.1
    pct_levels <- c("0.01", "0.05", "0.1")
    # Take subset of peaks based on fold change
    fc_subsets <- c(0.1, 0.3, 0.585, 1, 2)
    for(analysis_type in analysis_types) {
      print(analysis_type)
      analysis_type_dir <- paste0(cell_type_dir, analysis_type, "/")
      if (!dir.exists(analysis_type_dir)) {dir.create(analysis_type_dir, recursive = TRUE)}
      for(pct_level in pct_levels) {
        print(pct_level)
        pct_level_dir <- paste0(analysis_type_dir, pct_level, "/")
        if (!dir.exists(pct_level_dir)) {dir.create(pct_level_dir, recursive = TRUE)}
        no_bg_dir <- paste0(pct_level_dir, "without_bg", "/")
        if (!dir.exists(no_bg_dir)) {dir.create(no_bg_dir, recursive = TRUE)}
        bg_dir <- paste0(pct_level_dir, "with_bg", "/")
        if (!dir.exists(bg_dir)) {dir.create(bg_dir, recursive = TRUE)}
        for(fc_subset in fc_subsets) {
          print(fc_subset)
          # Read in peaks
          current_peaks <- read.table(paste0(motif_input_dir, "D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-time_point-controlling_for_subject_id_",
                                             analysis_type, "_pct_", pct_level, ".tsv"), 
                                      sep = "\t", header = TRUE)
          # Take subset of peaks based on fold change
          # Note that we filter FC based on single cell FC
          if(analysis_type == "sc") {
            pos_peaks <- current_peaks[current_peaks$avg_log2FC > 0,]
            neg_peaks <- current_peaks[current_peaks$avg_log2FC < 0,]
            pos_peaks <- pos_peaks[pos_peaks$avg_log2FC >= fc_subset,]
            neg_peaks <- neg_peaks[neg_peaks$avg_log2FC <= -fc_subset,]
          } else {
            pos_peaks <- current_peaks[current_peaks$sc_log2FC > 0,]
            neg_peaks <- current_peaks[current_peaks$sc_log2FC < 0,]
            pos_peaks <- pos_peaks[pos_peaks$sc_log2FC >= fc_subset,]
            neg_peaks <- neg_peaks[neg_peaks$sc_log2FC <= -fc_subset,]
          }
          # If we have no positive peaks, don't proceed with pos analysis
          if(nrow(pos_peaks) > 0) {
            total_pos_peaks <- nrow(pos_peaks)
            # Grab peak names 
            if(analysis_type == "sc") {
              pos_query_feature <- rownames(pos_peaks)
            } else {
              pos_query_feature <- pos_peaks$Peak_Name
            }
            # Find background peaks for pos and neg peaks
            pos.peaks.matched <- MatchRegionStats(
              meta.feature = meta.feature[open.peaks, ],
              query.feature = meta.feature[pos_query_feature, ],
              n = 40000
            )
            # Find motifs (pos, pos with bg, neg, neg with bg)
            pos_motifs <- FindMotifs(
              object = cells_subset,
              features = pos_query_feature
            )
            pos_motifs_with_bg <- FindMotifs(
              object = cells_subset,
              features = pos_query_feature,
              background = pos.peaks.matched
            )
            write.table(pos_motifs, file = paste0(no_bg_dir, "D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-",
                                                    analysis_type, "_pct_", pct_level, "_FC_", fc_subset, "_total_peaks_", total_pos_peaks, "_pos_motifs.tsv"), sep = "\t", quote = FALSE)
            write.table(pos_motifs_with_bg, file = paste0(bg_dir, "D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-",
                                                            analysis_type, "_pct_", pct_level, "_FC_", fc_subset, "_total_peaks_", total_pos_peaks, "_pos_motifs_with_bg.tsv"), sep = "\t", quote = FALSE)
          } 
          
          # If we have no negative peaks, don't proceed with neg analysis
          if(nrow(neg_peaks) > 0) {
            total_neg_peaks <- nrow(neg_peaks)
            # Grab peak names 
            if(analysis_type == "sc") {
              neg_query_feature <- rownames(neg_peaks)
            } else {
              neg_query_feature <- neg_peaks$Peak_Name
            }
            neg.peaks.matched <- MatchRegionStats(
              meta.feature = meta.feature[open.peaks, ],
              query.feature = meta.feature[neg_query_feature, ],
              n = 40000
            )
            neg_motifs <- FindMotifs(
              object = cells_subset,
              features = neg_query_feature
            )
            neg_motifs_with_bg <- FindMotifs(
              object = cells_subset,
              features = neg_query_feature,
              background = neg.peaks.matched
            )
            negative_FC_subset <- -fc_subset
            write.table(neg_motifs, file = paste0(no_bg_dir, "D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-",
                                                  analysis_type, "_pct_", pct_level, "_FC_", negative_FC_subset, "_total_peaks_", total_neg_peaks, "_neg_motifs.tsv"), sep = "\t", quote = FALSE)
            write.table(neg_motifs_with_bg, file = paste0(bg_dir, "D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-",
                                                          analysis_type, "_pct_", pct_level, "_FC_", negative_FC_subset, "_total_peaks_", total_neg_peaks, "_neg_motifs_with_bg.tsv"), sep = "\t", quote = FALSE)
          }
        }
      }
    }
  }
}
