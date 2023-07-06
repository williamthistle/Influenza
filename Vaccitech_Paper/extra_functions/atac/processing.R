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
  # Combine cell types
  Cell_type_combined = proj$predictedGroup
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
  proj$predictedGroup <- replace(proj$predictedGroup, proj$predictedGroup == "CD4 Naive", "T Naive")
  proj$predictedGroup <- replace(proj$predictedGroup, proj$predictedGroup == "CD8 Naive", "T Naive")
  proj$predictedGroup <- replace(proj$predictedGroup, proj$predictedGroup == "NK_CD56bright", "NK")
  proj$predictedGroup <- replace(proj$predictedGroup, proj$predictedGroup == "Treg", "T Naive")
  cM <- as.matrix(ArchR::confusionMatrix(proj$Harmony_clusters, proj$predictedGroup))
  Cell_type_voting <- proj$Harmony_clusters
  pre_cluster <- rownames(cM)
  max_celltype <- colnames(cM)[apply(cM, 1 , which.max)]
  for (m in c(1:length(pre_cluster))){
    idxSample <- which(proj$Harmony_clusters == pre_cluster[m])
    Cell_type_voting[idxSample] <- max_celltype[m]
  }
  proj <- addCellColData(ArchRProj = proj, data = Cell_type_voting, cells = proj$cellNames, name = "Cell_type_voting", force = TRUE)
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
  unique_cluster_ids <- unique(proj$Harmony_clusters)
  unique_cluster_ids <- unique_cluster_ids[order(nchar(unique_cluster_ids), unique_cluster_ids)]
  for (cluster in unique_cluster_ids) {
    idxPass <- which(proj$Harmony_clusters %in% cluster)
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
  idxPass <- which(proj$Harmony_clusters %in% cluster_identities)
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
}

# Method to create the cell type proportion file for MAGICAL
create_cell_type_proportion_MAGICAL_atac <- function(proj, analysis_dir, metadata_categories, day_metadata) {
  for(metadata_category in metadata_categories) {
    if(metadata_category == "day") {
      metadata <- day_metadata
    }
    cell_type_proportions_df <- data.frame("Condition" = metadata, "Sample_name" = paste0("Sample_", names(metadata)))
    total_cell_counts_df <- data.frame("Sample_name" = paste0("Sample_", names(metadata)))
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
        cell_type_proportions <- append(cell_type_proportions, cell_counts / total_cell_counts_df[total_cell_counts_df$Sample_name == paste0("Sample_", sample_id),]$cell_counts)
      }
      temp_df <- data.frame(cell_type_proportions)
      names(temp_df)[names(temp_df) == "cell_type_proportions"] <- cell_type
      cell_type_proportions_df <- cbind(cell_type_proportions_df, temp_df)
    }
    write.csv(cell_type_proportions_df, file = paste0(analysis_dir, "ATAC_cell_type_proportion_", metadata_category, ".csv"), quote = FALSE, row.names = FALSE)
  }
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
    cutOff= 0.05,
    force = TRUE
  )
  return(proj)
}

# Calculate differentially accessible peaks for each cell type
calculate_daps_for_each_cell_type <- function(proj, differential_peaks_dir) {
  # Calculate differential accessible peaks for each cell type
  for (cell_type in unique(proj$Cell_type_voting)) {
    print(cell_type)
    # Grab cells associated with cell type
    idxPass <- which(proj$Cell_type_voting %in% cell_type)
    print(length(idxPass))
    cellsPass <- proj$cellNames[idxPass]
    cells_subset <- proj[cellsPass,]
    print(cells_subset)
    # Find DAPs
    marker_D28_D1 <- getMarkerFeatures(ArchRProj = cells_subset, useMatrix = "PeakMatrix", groupBy = "day",
                                       testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"), normBy = "ReadsInPeaks", maxCells = 15000,
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
    write.xlsx(list_of_datasets, file = paste0(differential_peaks_dir, cell_type, "_", "D28_D1_diff.xlsx"), colNames = FALSE)
  }
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
}

# Create pseudobulk count files for each cell type
create_pseudobulk_atac <- function(proj, pseudobulk_dir) {
  # Create text file for each cell type containing pseudobulk counts for peaks
  Cell_types <- unique(final_proj$Cell_type_voting)
  sample.names <- unique(final_proj$Sample)
  peaks <- getPeakSet(proj)
  peak_count <- getMatrixFromProject(ArchRProj = proj, useMatrix = "PeakMatrix", useSeqnames = NULL, verbose = TRUE,binarize = FALSE,threads = getArchRThreads(),logFile = createLogFile("getMatrixFromProject"))
  for (i in c(1:length(Cell_types))){
    pseudo_bulk <- matrix(nrow = length(peaks), ncol = length(sample.names), 0)
    # We add Sample_ to the beginning of each column name to avoid MATLAB (MAGICAL) complaining
    colnames(pseudo_bulk)<-paste0("Sample_", sample.names)
    rownames(pseudo_bulk)<-peaks@elementMetadata$idx
    for (s in c(1:length(sample.names))){
      idxMatch <- which(str_detect(peak_count$Cell_type_voting,Cell_types[i]) & str_detect(as.character(peak_count$Sample),sample.names[s]))
      if (length(idxMatch)>1){
        pseudo_bulk[,s] = Matrix::rowSums(peak_count@assays@data$PeakMatrix[,idxMatch])
      }
    }
    cell_type <- sub(" ", "_", Cell_types[i])
    write.table(pseudo_bulk, file = paste0(pseudo_bulk_dir, "pseudo_bulk_ATAC_count_", cell_type, ".txt"), quote = FALSE, sep = "\t", col.names = NA)
  }
}