home_dir <- "~/SPEEDI"
source(paste0(home_dir, "/prototype_API.R"))

data_path <- "/data/home/wat2/high_vs_low_viral_load_D28/snRNA_seq_data/"
sample_id_list <- list.dirs(data_path, recursive = FALSE)
sample_id_list <- strsplit(sample_id_list, "/")
sample_id_list <- unlist(lapply(sample_id_list, tail, n = 1L))
output_dir <- "/data/home/wat2/high_vs_low_viral_load_D28/snRNA_seq_data_output/"

all_sc_exp_matrices <- Read_h5(data_path, sample_id_list)
sc_obj <- FilterRawData(all_sc_exp_matrices, human = TRUE)
rm(all_sc_exp_matrices)
sc_obj <- InitialProcessing(sc_obj, human = TRUE)
save.image(paste0(output_dir, "3_high_vs_low_viral_load_D28_placebo_multiome.RData"))
sc_obj <- InferBatches(sc_obj)
sc_obj <- IntegrateByBatch(sc_obj)
save.image(paste0(output_dir, "5_high_vs_low_viral_load_D28_placebo_multiome.RData"))
sc_obj <- VisualizeIntegration(sc_obj)
save.image(paste0(output_dir, "6_high_vs_low_viral_load_D28_placebo_multiome.RData"))
load(paste0(output_dir, "6_high_vs_low_viral_load_D28_placebo_multiome.RData"))
reference <- LoadReference("PBMC", human = TRUE)

### Test ###
# Remove certain cell types we're not interested in
#idx <- which(reference$celltype.l2 %in% c("Doublet", "B intermediate", "CD4 CTL", "gdT", "dnT", "ILC"))
#reference <- reference[,-idx]
#idx <- which(reference$celltype.l3 == "Treg Naive")
#reference <- reference[,-idx]
### Test ###


sc_obj <- MapCellTypes(sc_obj, reference)
# Print by cell type (majority vote)
DimPlot(sc_obj, reduction = "umap", group.by = "predicted_celltype_majority_vote", label = TRUE,
        label.size = 3, repel = TRUE, raster = FALSE) + 
  labs(title = "scRNA-seq Data Integration \n (X Samples, X Cells)") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(output_dir, paste0("high_vs_low_viral_load_D28_multiome_test", ".png")), device = "png", dpi = 300)
save.image(paste0(output_dir, "7_high_vs_low_viral_load_D28_placebo_multiome.RData"))

load(paste0(output_dir, "7_high_vs_low_viral_load_D28_placebo_multiome.RData"))
# Print by cluster number
DimPlot(sc_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
        label.size = 3, repel = TRUE, raster = FALSE) + 
  labs(title = "scRNA-seq Data Integration \n (X Samples, X Cells)") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(output_dir, paste0("high_vs_low_viral_load_D28_multiome_clusters", ".png")), device = "png", dpi = 300)

cell_names <- rownames(sc_obj@meta.data)
sc_obj <- AddMetaData(sc_obj, metadata = cell_names, col.name = "cell_name")

# Let's remove one cluster at a time and see how plot looks after removing each cluster. This will give us a better idea
# of which clusters are messy
for(cluster_id in unique(sc_obj$seurat_clusters)) {
  print(cluster_id)
  # Remove current cluster
  idxPass <- which(Idents(sc_obj) %in% c(cluster_id))
  cellsPass <- names(sc_obj$orig.ident[-idxPass])
  sc_obj.minus.clusters <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
  # Print cell type plot without cluster
  DimPlot(sc_obj.minus.clusters, reduction = "umap", group.by = "predicted_celltype_majority_vote", label = TRUE,
          label.size = 3, repel = TRUE, raster = FALSE) + 
    labs(title = "scRNA-seq Data Integration \n (X Samples, X Cells)") + 
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(output_dir, paste0("high_vs_low_viral_load_D28_multiome_without_cluster_", cluster_id, ".png")), device = "png", dpi = 300)
}

# Remove messy clusters - 10 is maybe questionable to remove (removes CD4 Naive)
idxPass <- which(Idents(sc_obj) %in% c(3, 10, 14))
cellsPass <- names(sc_obj$orig.ident[-idxPass])
sc_obj.minus.clusters <- subset(x = sc_obj, subset = cell_name %in% cellsPass)

# Reprint by cell type (majority vote)
DimPlot(sc_obj.minus.clusters, reduction = "umap", group.by = "predicted_celltype_majority_vote", label = TRUE,
        label.size = 3, repel = TRUE, raster = FALSE) + 
  labs(title = "scRNA-seq Data Integration \n (X Samples, X Cells)") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(output_dir, paste0("high_vs_low_viral_load_D28_multiome_without_messy_clusters.png")), device = "png", dpi = 300)

# Assign viral load (high or low) to each cell
high_viral_load <- c("b82bb7c75d47dac1", "3c4540710e55f7b1", "6f609a68dca1261f", "7b54cfac7e67b0fa")
low_viral_load <- c("abf6d19ee03be1e8", "216bb226181591dd", "d360f89cf9585dfe")
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

# Calculate cell type proportions 
cell_type_proportions_df <- data.frame("Condition" = viral_load_label, "Sample_name" = all_viral_load)
total_cell_counts_df <- data.frame("Condition" = viral_load_label, "Sample_name" = all_viral_load)
# Find total cell counts for each sample
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
write.csv(cell_type_proportions_df, file = paste0(output_dir, "RNA_cell_type_proportion.csv"), quote = FALSE, row.names = FALSE)
write.csv(total_cell_counts_df, file = paste0(output_dir, "RNA_cell_counts.csv"), quote = FALSE, row.names = FALSE)
#cell_type_proportions_plot <- ggplot(cell_type_proportions_df, aes(fill=condition, y=value, x=specie)) + 
#  geom_bar(position="stack", stat="identity")
  
#ggsave(paste0(output_dir, paste0("cell_type_proportions.png")), device = "png", dpi = 300)

# For one sample at a time
for(sample_id in sample_id_list) {
  all_sc_exp_matrices <- Read_h5(data_path, sample_id)
  sc_obj <- FilterRawData(all_sc_exp_matrices, human = TRUE)
  sc_obj <- InitialProcessing(sc_obj, human = TRUE)
  #sc_obj <- InferBatches(sc_obj)
  #sc_obj <- IntegrateByBatch(sc_obj)
  sc_obj <- VisualizeIntegration(sc_obj)
  reference <- LoadReference("PBMC", human = TRUE)
  sc_obj <- MapCellTypes(sc_obj, reference)
  DimPlot(sc_obj, reduction = "umap", group.by = "predicted_celltype_majority_vote", label = TRUE,
        label.size = 3, repel = TRUE, raster = FALSE) + 
  labs(title = "scRNA-seq Data Integration \n (X Samples, X Cells)") + 
  theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(output_dir, paste0(sample_id, ".png")), device = "png", dpi = 300)
}

#save.image(paste0(output_dir, "D28_placebo_multiome.RData"))
load(paste0(output_dir, "D28_placebo_multiome.RData"))