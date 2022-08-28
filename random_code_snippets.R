
output_dir <- "~/scRNA_seq_data_output/"
load(paste0(output_dir, "integrated_obj_after_predictions.RData"))

output_dir <- "~/scATAC_seq_data_output/"
load(paste0(output_dir, "atac_after_cell_type_voting.RData"))


for (res in seq(0, 1, 0.1)) {
  flu.combined.sct <- FindClusters(flu.combined.sct, resolution = res)
}

clustree(flu.combined.sct, prefix = "integrated_snn_res.")
ggsave(paste0(output_dir, "flu.combined.sct.cluster.trees.PDF"), device = "pdf")

proj <- saveArchRProject()
proj <- loadArchRProject(path = paste0(output_dir, "ArchR"))

load(paste0(output_dir, "atac_after_cell_type_voting_pre.RData"))

cluster_predictions <- vector()
cluster_distributions <- list()
idx <- 1
for (cluster in unique(proj.filtered$Clusters)) {
  idxPass <- which(proj.filtered$Clusters %in% cluster)
  cellsPass <- proj.filtered$cellNames[idxPass]
  filtered_cluster <-proj.filtered[cellsPass,]
  cluster_predictions <- append(cluster_predictions, table(filtered_cluster$Cell_type_voting))
  cluster_distributions[[idx]] <- table(filtered_cluster$Cell_type_combined)
  idx <- idx + 1
}

names(cluster_predictions) <- paste(unique(proj.filtered$Clusters), "-", names(cluster_predictions))


