dimensionality_reduc <- function(proj) {
  addArchRThreads(threads = 8)
  proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2,  force = TRUE,
                          clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 30),
                          varFeatures = 25000, dimsToUse = 2:30)
  proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", force = TRUE)
  proj <- addClusters(input = proj, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 5, knnAssign = 30, maxClusters = NULL, force = TRUE)
  return(proj)
}

plot_atac_after_filtering <- function(proj, date) {
  # UMAP plots colored by condition, sample, cluster ID, and TSS enrichment
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Conditions", embedding = "UMAP", force = TRUE)
  p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
  p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE)
  p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE)
  plotPDF(p1,p2,p3,p4, name = paste0("Integrated_Clustering_snRNA_", date), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
}

plot_atac_after_integration <- function(proj, date) {
  pal <- paletteDiscrete(values = proj$predictedGroup)
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP", pal = pal, force = TRUE)
  plotPDF(p1, name = paste0("Integrated_annotated_with_gene_integration_matrix_", date), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
}

plot_atac_after_majority_vote_or_subset <- function(proj, date) {
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE)
  p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
  p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", force = TRUE)
  p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE)
  plotPDF(p1,p2,p3,p4, name = paste0("Integrated_Clustering_Gene_Integration_Majority_Vote_or_Subset", date), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
}

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

add_rna_labels_for_atac_data <- function(proj, source_rna_file, use_rna_labels, subset_to_rna) {
  if(use_rna_labels) {
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
  }
  return(proj)
}

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
  return(proj)
}

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

get_cluster_info <- function(proj) {
  # See how clusters are distributed
  cluster_cell_type_predictions <- vector()
  cluster_cell_type_distributions <- list()
  cluster_sample_distributions <- list()
  cluster_viral_load_distributions <- list()
  cluster_day_distributions <- list()
  cluster_sex_distributions <- list()
  idx <- 1
  unique_cluster_ids <- unique(proj$Clusters)
  unique_cluster_ids <- unique_cluster_ids[order(nchar(unique_cluster_ids), unique_cluster_ids)]
  for (cluster in unique_cluster_ids) {
    idxPass <- which(proj$Clusters %in% cluster)
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

override_cluster_label <- function(proj, cluster_identities, cluster_label) {
  idxPass <- which(proj$Clusters %in% cluster_identities)
  proj$Cell_type_voting[idxPass] <- cluster_label
  return(proj)
}
