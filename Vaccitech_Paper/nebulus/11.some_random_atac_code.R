# Cluster info for ATAC clusters
get_cluster_info <- function(proj) {
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

# Majority vote for ATAC clusters
seurat_clusters <- as.factor(atac_proj$Clusters)
predictedGroup <- atac_proj$predictedGroup
predictedScore <- atac_proj$predictedScore
votes <- c()
vote_levels <- as.character(levels(seurat_clusters))
for (i in vote_levels) {
  cluster_cells <- which(atac_proj$Clusters == i)
  sub_predictedGroup <- predictedGroup[cluster_cells]
  sub_predictedScore <- predictedScore[cluster_cells]
  gmeans <- c()
  cell_types <- c()
  for (j in names(prop.table(table(sub_predictedGroup))[prop.table(table(sub_predictedGroup)) > 
                                                        0.25])) {
    cell_types <- c(cell_types, j)
    cell_type_cells <- which(sub_predictedGroup == j)
    gmeans <- c(gmeans, exp(mean(log(sub_predictedScore[cell_type_cells]))))
  }
  if (!is.null(cell_types)) {
    vote_levels[vote_levels %in% as.character(i)] <- cell_types[which.max(gmeans)]
  }
  else {
    vote_levels[vote_levels %in% as.character(i)] <- "Undetermined"
  }
}
cell_type_voting <- seurat_clusters
levels(cell_type_voting) <- vote_levels
cell_type_voting <- as.character(cell_type_voting)
atac_proj <- ArchR::addCellColData(ArchRProj = atac_proj, data = cell_type_voting, 
                              cells = atac_proj$cellNames, name = "Cell_type_voting", force = TRUE)

# Override ATAC cluster
override_cluster_label_atac <- function(proj, cluster_identities, cluster_label) {
  idxPass <- which(proj$Clusters %in% cluster_identities)
  proj$Cell_type_voting[idxPass] <- cluster_label
  return(proj)
}


p1 <- ArchR::plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", 
                           name = "Clusters", embedding = "UMAP", 
                           force = TRUE, keepAxis = TRUE)
ggplot2::ggsave(filename = paste0(output_dir, "Xi_By_Cluster.png"), 
                plot = p1, device = "png", width = 8, height = 8, 
                units = "in")








# Subsetting clusters 
idxPass <- which(atac_proj$Clusters %in% c("C1", "C2", "C3", "C8", "C11", "C12", "C13", "C14", "C24", "C28", "C36"))
cellsPass <- atac_proj$cellNames[-idxPass]
atac_proj_minus_clusters <- atac_proj[cellsPass, ]

atac_proj_minus_clusters <- override_cluster_label_atac(atac_proj_minus_clusters, c("C16"), "CD8 Memory")
atac_proj_minus_clusters <- override_cluster_label_atac(atac_proj_minus_clusters, c("C20"), "Proliferating")
atac_proj_minus_clusters <- override_cluster_label_atac(atac_proj_minus_clusters, c("C25", "C34"), "CD4 Naive")


p1 <- ArchR::plotEmbedding(ArchRProj = atac_proj_minus_clusters, colorBy = "cellColData", 
                           name = "Cell_type_voting", embedding = "UMAP", 
                           force = TRUE, keepAxis = TRUE)

ggplot2::ggsave(filename = paste0(output_dir, "Xi_cell_type_voting_minus_clusters_new.png"), 
                plot = p1, device = "png", width = 8, height = 8, 
                units = "in")

p1 <- ArchR::plotEmbedding(ArchRProj = atac_proj_minus_clusters, colorBy = "cellColData", 
                           name = "Clusters", embedding = "UMAP", 
                           force = TRUE, keepAxis = TRUE)

ggplot2::ggsave(filename = paste0(output_dir, "Xi_cell_type_voting_minus_clusters_clusters.png"), 
                plot = p1, device = "png", width = 8, height = 8, 
                units = "in")