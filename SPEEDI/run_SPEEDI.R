library(clustree)

home_dir <- "~/SPEEDI"
source(paste0(home_dir, "/prototype_API.R"))

# naming_token is used to name output files
naming_token <- "all_high_vs_low_viral_load_D28_reduced_references"

# data_path is where input data are stored
data_path <- paste0("/data/home/wat2/", naming_token, "/RNA_seq_data/")
# Get list of samples that will be processed
sample_id_list <- list.dirs(data_path, recursive = FALSE)
sample_id_list <- strsplit(sample_id_list, "/")
sample_id_list <- unlist(lapply(sample_id_list, tail, n = 1L))
output_dir <- paste0("/data/home/wat2/", naming_token, "/RNA_seq_data_output/")

# Run SPEEDI on our list of samples
save_progress <- TRUE
sc_obj <- run_SPEEDI(data_path, output_dir, sample_id_list, naming_token, save_progress, use_simplified_reference = TRUE)
if(!save_progress) {
  save.image(paste0(output_dir, "7_", naming_token, ".RData"))
}

#load(paste0(output_dir, "7_", naming_token, ".RData"))

# Print UMAP by cell type (majority vote) and by cluster number - it will currently be messy
sample_count <- length(sample_id_list)
cell_count <- length(sc_obj$cell_name)
current_title <- paste0("scRNA-seq and/or snRNA-seq Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
print_UMAP(sc_obj, "predicted_celltype_majority_vote", current_title, output_dir, naming_token, "_clusters_by_cell_type_messy.png")
print_UMAP(sc_obj, "seurat_clusters", current_title, output_dir, naming_token, "_clusters_by_cluster_num_messy.png")

# Let's use clustree to try to figure out the best clustering resolution
for (res in seq(3, 6, 0.3)) {
  sc_obj <- FindClusters(sc_obj, resolution = res)
}
clustree(sc_obj, prefix = "integrated_snn_res.3")
ggsave(paste0(output_dir, naming_token, "_cluster.trees.PNG"), device = "png", width = 8, height = 8, units = "in")

# Now, we should re-run our majority vote with the correct resolution
best_res <- 3
associated_res_attribute <- paste0("integrated_snn_res.", best_res)
sc_obj <- MajorityVote(sc_obj, best_res, associated_res_attribute)

# To decide which clusters we need to remove, we will use two tactics:
# 1. Capture information about cell type distributions in each cluster (and cluster mean S score, G2M score, and CC diff for cell cycle info)
#    High S score and high G2M score seem to indicate Proliferating cluster
# 2. Remove one cluster at a time and see how plot looks after removing each cluster
# 1
cluster_distributions <- list()
cluster_predictions <- vector()
cluster_mean_S_score <- vector()
cluster_mean_G2M_score <- vector()
cluster_mean_CC_difference <- vector()
cluster_ids <- vector()
idx <- 1
for (cluster in levels(sc_obj)) {
  idxPass <- which(Idents(sc_obj) %in% cluster)
  cellsPass <- names(sc_obj$orig.ident[idxPass])
  filtered_cluster <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
  cluster_prediction <- sort(table(filtered_cluster$predicted_celltype_majority_vote), decreasing = TRUE)[1]
  cluster_predictions <- append(cluster_predictions, cluster_prediction)
  cluster_distributions[[idx]] <- table(filtered_cluster$predicted.id)
  cluster_mean_S_score <- append(cluster_mean_S_score, mean(filtered_cluster$S.Score))
  cluster_mean_G2M_score <- append(cluster_mean_G2M_score, mean(filtered_cluster$G2M.Score))
  cluster_mean_CC_difference <- append(cluster_mean_CC_difference, mean(filtered_cluster$CC.Difference))
  cluster_ids <- append(cluster_ids, cluster)
  idx <- idx + 1
}
names(cluster_predictions) <- paste(levels(sc_obj), "-", names(cluster_predictions))
cell_cycle_df <- data.frame("Cluster" = cluster_ids, "S" = cluster_mean_S_score, "G2M" = cluster_mean_G2M_score, "CC Diff" = cluster_mean_CC_difference) 
#2
for(cluster_id in unique(sc_obj$seurat_clusters)) {
  print(cluster_id)
  # Remove current cluster
  idxPass <- which(Idents(sc_obj) %in% c(cluster_id))
  cellsPass <- names(sc_obj$orig.ident[-idxPass])
  sc_obj.minus.current.cluster <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
  cell_count <- length(sc_obj.minus.current.cluster$cell_name)
  current_title <- paste0("scRNA-seq and/or snRNA-seq Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
  # Print cell type plot without cluster
  print_UMAP(sc_obj.minus.current.cluster, "predicted_celltype_majority_vote", current_title, output_dir, naming_token, paste0("_clusters_by_cell_type_without_cluster_", cluster_id, ".png"))
}

# Remove messy clusters
messy_clusters <- c(1, 3, 8, 16, 19, 20, 29)
idxPass <- which(Idents(sc_obj) %in% messy_clusters)
cellsPass <- names(sc_obj$orig.ident[-idxPass])
sc_obj.minus.messy.clusters <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
# Print cell type plot without messy clusters
cell_count <- length(sc_obj.minus.current.cluster$cell_name)
current_title <- paste0("scRNA-seq and/or snRNA-seq Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
print_UMAP(sc_obj.minus.messy.clusters, "predicted_celltype_majority_vote", current_title, output_dir, naming_token, "_clusters_by_cell_type_without_messy_clusters.png")

# Assign viral load (high or low) to each cell
high_viral_load <- c("344465e6c3a8cb53", "4534496c580cb408", "3731a6247ae23831", "b82bb7c75d47dac1", "3c4540710e55f7b1", "6f609a68dca1261f", "7b54cfac7e67b0fa")
low_viral_load <- c("a7e1ca56bfd75e11", "abf6d19ee03be1e8", "216bb226181591dd", "d360f89cf9585dfe")
all_viral_load <- c(high_viral_load, low_viral_load)
viral_load_label <- c(rep("HIGH", length(high_viral_load)), rep("LOW", length(low_viral_load)))
viral_load_vec <- c()
for(current_sample in sc_obj.minus.clusters$sample) {
  if(current_sample %in% high_viral_load) {
    viral_load_vec <- c(viral_load_vec, "HIGH")
  } else {
    viral_load_vec <- c(viral_load_vec, "LOW")
  }
}
sc_obj.minus.clusters$viral_load <- viral_load_vec

# Print viral load plot
print_UMAP(sc_obj.minus.messy.clusters, "viral_load", current_title, output_dir, naming_token, "_clusters_by_viral_load_without_messy_clusters.png")

# Calculate cell type proportions and cell counts for each cell type (split by sample)
cell_type_proportions_df <- data.frame("Condition" = viral_load_label, "Sample_name" = all_viral_load)
total_cell_counts_df <- data.frame("Condition" = viral_load_label, "Sample_name" = all_viral_load)
total_cell_counts <- c()
for (sample_id in all_viral_load) {
  idxPass <- which(sc_obj.minus.clusters$sample %in% sample_id)
  cellsPass <- names(sc_obj.minus.clusters$orig.ident[idxPass])
  sample_subset <- subset(x = sc_obj.minus.clusters, subset = cell_name %in% cellsPass)
  total_cell_counts <- c(total_cell_counts, ncol(sample_subset))
}
total_cell_counts_df <- cbind(total_cell_counts_df, total_cell_counts)

for (cell_type in unique(sc_obj.minus.clusters$predicted_celltype_majority_vote)) {
  cell_type_proportions <- vector()
  cell_counts <- c()
  print(cell_type)
  # Grab cells associated with cell type
  idxPass <- which(sc_obj.minus.clusters$predicted_celltype_majority_vote %in% cell_type)
  print(length(idxPass))
  cellsPass <- names(sc_obj.minus.clusters$orig.ident[idxPass])
  cells_subset <- subset(x = sc_obj.minus.clusters, subset = cell_name %in% cellsPass)
  for (sample_id in all_viral_load) {
    # Subset further based on cells associated with sample ID
    idxPass <- which(cells_subset$sample %in% sample_id)
    cellsPass <- names(cells_subset$orig.ident[idxPass])
    if (length(cellsPass) == 0) {
      cell_type_proportions <- append(cell_type_proportions, 0)
    } else {
      sample_subset <- subset(x = cells_subset, subset = cell_name %in% cellsPass)
      cell_count <- ncol(sample_subset)
      cell_counts <- c(cell_counts, cell_count)
      cell_type_proportions <- append(cell_type_proportions, cell_count / total_cell_counts_df[total_cell_counts_df$Sample_name == sample_id,]$total_cell_counts)
    }
  }
  temp_df <- data.frame(cell_counts)
  names(temp_df)[names(temp_df) == "cell_counts"] <- cell_type
  total_cell_counts_df <- cbind(total_cell_counts_df, temp_df)
  temp_df <- data.frame(cell_type_proportions)
  names(temp_df)[names(temp_df) == "cell_type_proportions"] <- cell_type
  cell_type_proportions_df <- cbind(cell_type_proportions_df, temp_df)
}
write.csv(cell_type_proportions_df, file = paste0(output_dir, naming_token, "_RNA_cell_type_proportion.reduced.references.csv"), quote = FALSE, row.names = FALSE)
write.csv(total_cell_counts_df, file = paste0(output_dir, naming_token, "_RNA_cell_counts.reduced.references.csv"), quote = FALSE, row.names = FALSE)

