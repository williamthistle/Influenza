library(clustree)

home_dir <- "~/SPEEDI"
source(paste0(home_dir, "/prototype_API.R"))

# naming_token is used to name output files
naming_token <- "high_vs_low_viral_load_D28"
#naming_token <- "all_high_vs_low_viral_load_D28_reduced_references"

# data_path is where input data are stored
data_path <- paste0("/data/home/wat2/", naming_token, "/RNA_seq_data/")
# Get list of samples that will be processed
sample_id_list <- list.dirs(data_path, recursive = FALSE)
sample_id_list <- strsplit(sample_id_list, "/")
sample_id_list <- unlist(lapply(sample_id_list, tail, n = 1L))
sample_count <- length(sample_id_list)
output_dir <- paste0("/data/home/wat2/", naming_token, "/RNA_seq_data_output/")

# Run SPEEDI on our list of samples
save_progress <- TRUE
sc_obj <- run_SPEEDI(data_path, output_dir, sample_id_list, naming_token, save_progress, use_simplified_reference = TRUE)
if(!save_progress) {
  save.image(paste0(output_dir, "7_", naming_token, ".RData"))
}

#load(paste0(output_dir, "7_", naming_token, ".RData"))

# Add cell names to Seurat object
cell_names <- rownames(sc_obj@meta.data)
sc_obj <- AddMetaData(sc_obj, metadata = cell_names, col.name = "cell_name")
# Assign viral load (high or low) to each cell
high_viral_load <- c("b82bb7c75d47dac1", "3c4540710e55f7b1", "6f609a68dca1261f", "7b54cfac7e67b0fa")
low_viral_load <- c("abf6d19ee03be1e8", "216bb226181591dd", "d360f89cf9585dfe")
all_viral_load <- c(high_viral_load, low_viral_load)
viral_load_label <- c(rep("HIGH", length(high_viral_load)), rep("LOW", length(low_viral_load)))
viral_load_vec <- c()
for(current_sample in sc_obj$sample) {
  if(current_sample %in% high_viral_load) {
    viral_load_vec <- c(viral_load_vec, "HIGH")
  } else {
    viral_load_vec <- c(viral_load_vec, "LOW")
  }
}
sc_obj$viral_load <- viral_load_vec

# Print UMAP by cell type (majority vote) and by cluster number - it will currently be messy
print_UMAP(sc_obj, sample_count, "predicted_celltype_majority_vote", output_dir, naming_token, "_clusters_by_cell_type-2-01.png")
print_UMAP(sc_obj, sample_count, "seurat_clusters", output_dir, naming_token, "_clusters_by_cluster_num-2-01.png")

# Print viral load plot
print_UMAP(sc_obj, sample_count, "viral_load", output_dir, naming_token, "_clusters_by_viral_load-1-31.png")

# Let's see what the cell type proportions are for our raw data (using predicted.id)
# We aren't that concerned about what the majority vote says right now because we still need to clean up our clusters
raw_cell_tables <- calculate_props_and_counts(sc_obj, viral_load_label, all_viral_load)
write.csv(raw_cell_tables[1], file = paste0(output_dir, naming_token, "_raw_RNA_cell_type_proportion_1-31.csv"), quote = FALSE, row.names = FALSE)
write.csv(raw_cell_tables[2], file = paste0(output_dir, naming_token, "_raw_RNA_cell_counts_1-31.csv"), quote = FALSE, row.names = FALSE)

# Let's also use Vincy's code to create a nice heatmap of cell type proportions
sc_obj$sample <- factor(sc_obj$sample, levels = c("3c4540710e55f7b1", "6f609a68dca1261f", "7b54cfac7e67b0fa", "b82bb7c75d47dac1", "216bb226181591dd", "abf6d19ee03be1e8", "d360f89cf9585dfe"))
raw_cell_type_proportion  <- as.matrix(table(sc_obj$sample, sc_obj$predicted.id))
raw_cell_type_proportion <- apply(raw_cell_type_proportion, 1, function(x){x/sum(x)})
# Combine B Cells
rownames(raw_cell_type_proportion)[2] <- "B Cells"
raw_cell_type_proportion[2,] <- raw_cell_type_proportion[2,] + raw_cell_type_proportion[3,] + raw_cell_type_proportion[4,]
raw_cell_type_proportion <- raw_cell_type_proportion[-c(3, 4), ]
# Combine Monocytes
rownames(raw_cell_type_proportion)[3] <- "Monocytes"
raw_cell_type_proportion[3,] <- raw_cell_type_proportion[3,] + raw_cell_type_proportion[4,]
raw_cell_type_proportion <- raw_cell_type_proportion[-c(4), ]
# Combine T Cells
rownames(raw_cell_type_proportion)[4] <- "T Cells"
raw_cell_type_proportion[4,] <- raw_cell_type_proportion[4,] + raw_cell_type_proportion[6,] + raw_cell_type_proportion[7,] + raw_cell_type_proportion[8,] + raw_cell_type_proportion[10,]
raw_cell_type_proportion <- raw_cell_type_proportion[-c(6, 7, 8, 10), ]
# Combine Rare T Cells
rownames(raw_cell_type_proportion)[5] <- "Rare T Cells"
raw_cell_type_proportion[5,] <- raw_cell_type_proportion[5,] + raw_cell_type_proportion[6,] + raw_cell_type_proportion[8,] + raw_cell_type_proportion[10,] + raw_cell_type_proportion[19,]
raw_cell_type_proportion <- raw_cell_type_proportion[-c(6, 8, 10, 19), ]
# Combine NK / ILC
rownames(raw_cell_type_proportion)[9] <- "NK / ILC"
raw_cell_type_proportion[9,] <- raw_cell_type_proportion[9,] + raw_cell_type_proportion[11,] + raw_cell_type_proportion[12,]
raw_cell_type_proportion <- raw_cell_type_proportion[-c(11, 12), ]

my_sample_col <- data.frame(sample = rep(c("HIGH", "LOW"), c(4,3))) #The study group by their order
row.names(my_sample_col) <- colnames(raw_cell_type_proportion)

library(pheatmap)
output.plot <- pheatmap(raw_cell_type_proportion, annotation_col = my_sample_col, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.3f")
ggsave(paste0(output_dir, naming_token, "_raw_cell_type_proportions_heatmap.PNG"), plot = output.plot, device = "png", width = 8, height = 8, units = "in")

# Let's use clustree to try to figure out the best clustering resolution
for (res in seq(0, 3, 0.3)) {
  sc_obj <- FindClusters(sc_obj, resolution = res)
}
clustree(sc_obj, prefix = "integrated_snn_res.")
ggsave(paste0(output_dir, naming_token, "_cluster.trees-1-31.PNG"), device = "png", width = 8, height = 8, units = "in")

# Now, we should re-run our majority vote with the correct resolution
best_res <- 1.5
sc_obj <- MajorityVote(sc_obj, best_res)

# To decide which clusters we need to remove, we will use two tactics:
# 1. Capture information about cell type distributions in each cluster (and cluster mean S score, G2M score, and CC diff for cell cycle info)
#    High S score and high G2M score seem to indicate Proliferating cluster
# 2. Remove one cluster at a time and see how plot looks after removing each cluster
# 1
raw_cluster_info <- capture_cluster_info(sc_obj)
# 2
for(cluster_id in unique(sc_obj$seurat_clusters)) {
  print(cluster_id)
  # Remove current cluster
  idxPass <- which(Idents(sc_obj) %in% c(cluster_id))
  cellsPass <- names(sc_obj$orig.ident[-idxPass])
  sc_obj.minus.current.cluster <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
  cell_count <- length(sc_obj.minus.current.cluster$cell_name)
  current_title <- paste0("scRNA-seq and/or snRNA-seq Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
  # Print cell type plot without cluster
  print_UMAP(sc_obj.minus.current.cluster, sample_count, "predicted_celltype_majority_vote", current_title, output_dir, naming_token, paste0("_clusters_by_cell_type_without_cluster_", cluster_id, "_02-01.png"))
}

# To confirm the cell type associated with a given cluster, we can also find the cluster markers (example below)
cluster43.markers <- FindMarkers(sc_obj, ident.1 = 43, min.pct = 0.25)

# Remove messy clusters
messy_clusters <- c(1, 3, 8, 9, 16, 19, 20, 28) # For multiome F
#messy_clusters <- c(9, 29, 40, 45, 47, 49, 50, 58, 60) # For multiome + scRNA-seq
idxPass <- which(Idents(sc_obj) %in% messy_clusters)
cellsPass <- names(sc_obj$orig.ident[-idxPass])
sc_obj.minus.messy.clusters <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
# The ideal resolution may have changed. Let's check now!
for (res in seq(0, 3, 0.3)) {
  sc_obj.minus.messy.clusters <- FindClusters(sc_obj.minus.messy.clusters, resolution = res)
}
clustree(sc_obj.minus.messy.clusters, prefix = "integrated_snn_res.")
ggsave(paste0(output_dir, naming_token, "_cluster.trees.minus.messy.clusters.PNG"), device = "png", width = 8, height = 8, units = "in")

# Redo majority vote with updated resolution - here, 1.2 looks the best
sc_obj.minus.messy.clusters <- MajorityVote(sc_obj.minus.messy.clusters)
# Print cell type plot without messy clusters
cell_count <- length(sc_obj.minus.messy.clusters$cell_name)
current_title <- paste0("scRNA-seq and/or snRNA-seq Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
print_UMAP(sc_obj.minus.messy.clusters, sample_count, "predicted_celltype_majority_vote", current_title, output_dir, naming_token, "_clusters_by_cell_type_without_messy_clusters.png")

# Update cluster predictions and produce new plots
processed_cluster_info <- capture_cluster_info(sc_obj.minus.messy.clusters)
#2
for(cluster_id in unique(sc_obj.minus.messy.clusters$seurat_clusters)) {
  print(cluster_id)
  # Remove current cluster
  idxPass <- which(Idents(sc_obj.minus.messy.clusters) %in% c(cluster_id))
  cellsPass <- names(sc_obj.minus.messy.clusters$orig.ident[-idxPass])
  sc_obj.minus.current.cluster <- subset(x = sc_obj.minus.messy.clusters, subset = cell_name %in% cellsPass)
  cell_count <- length(sc_obj.minus.current.cluster$cell_name)
  current_title <- paste0("scRNA-seq and/or snRNA-seq Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
  # Print cell type plot without cluster
  print_UMAP(sc_obj.minus.current.cluster, sample_count, "predicted_celltype_majority_vote", current_title, output_dir, naming_token, paste0("_clusters_by_cell_type_minus_messy_clusters_without_cluster_", cluster_id, ".png"))
}

# Remove some more messy clusters, I guess
messy_clusters_round_2 <- c(1, 11)
idxPass <- which(Idents(sc_obj.minus.messy.clusters) %in% messy_clusters_round_2)
cellsPass <- names(sc_obj.minus.messy.clusters$orig.ident[-idxPass])
sc_obj.minus.messy.clusters <- subset(x = sc_obj.minus.messy.clusters, subset = cell_name %in% cellsPass)

# Plot again
cell_count <- length(sc_obj.minus.messy.clusters$cell_name)
current_title <- paste0("scRNA-seq and/or snRNA-seq Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
print_UMAP(sc_obj.minus.messy.clusters, sample_count, "predicted_celltype_majority_vote", current_title, output_dir, naming_token, "_clusters_by_cell_type_without_messy_clusters_round_2.png")

# Print viral load plot
print_UMAP(sc_obj.minus.messy.clusters, "viral_load", current_title, output_dir, naming_token, "_clusters_by_viral_load_without_messy_clusters.png")

# Calculate cell type proportions and cell counts for each cell type (split by sample)
processed_cell_tables <- calculate_props_and_counts(sc_obj.minus.messy.clusters, viral_load_label, all_viral_load)
write.csv(processed_cell_tables[1], file = paste0(output_dir, naming_token, "_RNA_cell_type_proportion.csv"), quote = FALSE, row.names = FALSE)
write.csv(processed_cell_tables[2], file = paste0(output_dir, naming_token, "_RNA_cell_counts.csv"), quote = FALSE, row.names = FALSE)





