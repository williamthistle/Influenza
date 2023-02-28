for (i in 1:10) {
  SEED <- i
  set.seed(SEED)
  all_sc_exp_matrices <- Read_h5(data_path, sample_id_list)
  assign("sc_obj", FilterRawData(all_sc_exp_matrices, human = TRUE, remove_doublets = TRUE), envir = .GlobalEnv)
  rm(all_sc_exp_matrices)
  assign("sc_obj", InitialProcessing(sc_obj, human = TRUE), envir = .GlobalEnv)
  assign("sc_obj", InferBatches(sc_obj), envir = .GlobalEnv)
  assign("sc_obj", IntegrateByBatch(sc_obj), envir = .GlobalEnv)
  assign("sc_obj", VisualizeIntegration(sc_obj), envir = .GlobalEnv)
  reference <- LoadReference("PBMC", human = TRUE)
  assign("sc_obj", MapCellTypes(sc_obj, reference), envir = .GlobalEnv)
  rm(reference)
  cell_names <- rownames(sc_obj@meta.data)
  assign("sc_obj", AddMetaData(sc_obj, metadata = cell_names, col.name = "cell_name"), envir = .GlobalEnv)
  print_UMAP(sc_obj, sample_count, "predicted_celltype_majority_vote", output_dir, naming_token, paste0("_clusters_by_cell_type_majority_vote_", date, SEED, "_SPEEDI.png"))
  print_UMAP(sc_obj, sample_count, "seurat_clusters", output_dir, naming_token, paste0("_clusters_by_cluster_num_", date, SEED, "_SPEEDI.png"))
  print_UMAP(sc_obj, sample_count, "predicted.id", output_dir, naming_token, paste0("_clusters_by_cell_type_", date, SEED, "_SPEEDI.png"))
  save(sc_obj, file = paste0(output_dir, "7_", naming_token, "_sc_obj_with_doublets_intact_SPEEDI.rds"))
}

# Maybe I should process all the cells and then look at the doublet scores for the ones in the bridge clusters?
SEED <- 3
set.seed(SEED)
all_sc_exp_matrices <- Read_h5(data_path, sample_id_list)
assign("sc_obj_3", FilterRawData(all_sc_exp_matrices, human = TRUE, remove_doublets = TRUE), envir = .GlobalEnv)
rm(all_sc_exp_matrices)
SEED <- 4
set.seed(SEED)
all_sc_exp_matrices <- Read_h5(data_path, sample_id_list)
assign("sc_obj_4", FilterRawData(all_sc_exp_matrices, human = TRUE, remove_doublets = TRUE), envir = .GlobalEnv)
rm(all_sc_exp_matrices)
cell_names <- rownames(sc_obj_3@meta.data)
assign("sc_obj_3", AddMetaData(sc_obj_3, metadata = cell_names, col.name = "cell_name"), envir = .GlobalEnv)
cell_names <- rownames(sc_obj_4@meta.data)
assign("sc_obj_4", AddMetaData(sc_obj_4, metadata = cell_names, col.name = "cell_name"), envir = .GlobalEnv)

doublet_sc_obj_3 <- subset(x = sc_obj_3, subset = scDblFinder.class %in% "doublet")
doublet_sc_obj_4 <- subset(x = sc_obj_4, subset = scDblFinder.class %in% "doublet")
non_overlapping_cells_in_sc_obj_3 <- setdiff(doublet_sc_obj_3$cell_name, doublet_sc_obj_4$cell_name)
non_overlapping_cells_in_sc_obj_4 <- setdiff(doublet_sc_obj_4$cell_name, doublet_sc_obj_3$cell_name)


# How did SEED 3 score cells that were found to be doublets in SEED 3 but not SEED 4?
sc_obj_3_doublets_only_3 <- subset(x = sc_obj_3, subset = cell_name %in% non_overlapping_cells_in_sc_obj_3)
mean(sc_obj_3_doublets_only_3$scDblFinder.score)
median(sc_obj_3_doublets_only_3$scDblFinder.score)
min(sc_obj_3_doublets_only_3$scDblFinder.score)
max(sc_obj_3_doublets_only_3$scDblFinder.score)
# How did SEED 4 score cells that were found to be doublets in SEED 3 but not SEED 4?
sc_obj_4_doublets_only_3 <- subset(x = sc_obj_4, subset = cell_name %in% non_overlapping_cells_in_sc_obj_3)
mean(sc_obj_4_doublets_only_3$scDblFinder.score)
median(sc_obj_4_doublets_only_3$scDblFinder.score)
min(sc_obj_4_doublets_only_3$scDblFinder.score)
max(sc_obj_4_doublets_only_3$scDblFinder.score)
# How did SEED 3 score cells that were found to be doublets in SEED 4 but not SEED 3?
sc_obj_3_doublets_only_4 <- subset(x = sc_obj_3, subset = cell_name %in% non_overlapping_cells_in_sc_obj_4)
mean(sc_obj_3_doublets_only_4$scDblFinder.score)
median(sc_obj_3_doublets_only_4$scDblFinder.score)
min(sc_obj_3_doublets_only_4$scDblFinder.score)
max(sc_obj_3_doublets_only_4$scDblFinder.score)
# How did SEED 4 score cells that were found to be doublets in SEED 4 but not SEED 3?
sc_obj_4_doublets_only_4 <- subset(x = sc_obj_4, subset = cell_name %in% non_overlapping_cells_in_sc_obj_4)
mean(sc_obj_4_doublets_only_4$scDblFinder.score)
median(sc_obj_4_doublets_only_4$scDblFinder.score)
min(sc_obj_4_doublets_only_4$scDblFinder.score)
max(sc_obj_4_doublets_only_4$scDblFinder.score)



# SC3
#Step 2: Filtering out bad samples...
#Removing doublets...
#Number of doublets removed in each sample:

#216bb226181591dd 3c4540710e55f7b1 6f609a68dca1261f 7b54cfac7e67b0fa
#2815             3683             3998             4370
#abf6d19ee03be1e8 b82bb7c75d47dac1 d360f89cf9585dfe
#2381             1575             3832

# SC4
#Number of doublets removed in each sample:

#216bb226181591dd 3c4540710e55f7b1 6f609a68dca1261f 7b54cfac7e67b0fa
#2733             3803             4031             4379
#abf6d19ee03be1e8 b82bb7c75d47dac1 d360f89cf9585dfe
#2249             1603             3776

# > length(sc_obj_3$cell_name)
# [1] 53018
# > length(sc_obj_4$cell_name)
# [1] 53032
# > length(intersect(sc_obj_3$cell_name, sc_obj_4$cell_name))
# [1] 51272


# Calculate what it looks like after removing increasing numbers of doublets
for (i in seq(0.9, 0.1, by=-0.1)) {
  sc_obj_doublets_removed <- subset(x = sc_obj, subset = scDblFinder.score < i)
  sc_obj_doublets_removed <- ScaleData(sc_obj_doublets_removed, verbose = T)
  sc_obj_doublets_removed <- RunPCA(sc_obj_doublets_removed, npcs = 30, approx = T, verbose = T, seed.use = SEED)
  sc_obj_doublets_removed <- RunUMAP(sc_obj_doublets_removed, reduction = "pca", dims = 1:30, seed.use = SEED, return.model = T)
  sc_obj_doublets_removed <- MapCellTypes(sc_obj_doublets_removed, reference)
  print_UMAP(sc_obj_doublets_removed, sample_count, "predicted_celltype_majority_vote", output_dir, naming_token, paste0("_clusters_by_cell_type_majority_vote_", date, "_", i, ".png"))
  print_UMAP(sc_obj_doublets_removed, sample_count, "seurat_clusters", output_dir, naming_token, paste0("_clusters_by_cluster_num_", date, "_", i, ".png"))
  print_UMAP(sc_obj_doublets_removed, sample_count, "predicted.id", output_dir, naming_token, paste0("_clusters_by_cell_type_", date, "_", i, ".png"))
}

# See overlap between ATAC and RNA-seq
filtered_ATAC_cells <- read.table(paste0(home_dir, "/ATAC_filtered_cells.txt"), comment.char = "")
filtered_ATAC_cells <- filtered_ATAC_cells$V1
sc_obj_overlap_ATAC <- subset(sc_obj, cell_name %in% filtered_ATAC_cells)
# 36129 remain after overlap
sc_obj_overlap_ATAC <- RunPCA(sc_obj_overlap_ATAC, npcs = 30, approx = T, verbose = T, seed.use = SEED)
sc_obj_overlap_ATAC <- RunUMAP(sc_obj_overlap_ATAC, reduction = "pca", dims = 1:30, seed.use = SEED, return.model = T)
sc_obj_overlap_ATAC <- MajorityVote(sc_obj_overlap_ATAC)
print_UMAP(sc_obj_overlap_ATAC, sample_count, "predicted_celltype_majority_vote", output_dir, naming_token, paste0("_clusters_by_cell_type_majority_vote_", date, "_ATAC_overlap.png"))
for (i in seq(0.9, 0.1, by=-0.1)) {
  sc_obj_overlap_ATAC_doublets_removed <- subset(x = sc_obj_overlap_ATAC, subset = scDblFinder.score < i)
  sc_obj_overlap_ATAC_doublets_removed <- RunPCA(sc_obj_overlap_ATAC_doublets_removed, npcs = 30, approx = T, verbose = T, seed.use = SEED)
  sc_obj_overlap_ATAC_doublets_removed <- RunUMAP(sc_obj_overlap_ATAC_doublets_removed, reduction = "pca", dims = 1:30, seed.use = SEED, return.model = T)
  sc_obj_overlap_ATAC_doublets_removed <- MajorityVote(sc_obj_overlap_ATAC_doublets_removed)
  print_UMAP(sc_obj_overlap_ATAC_doublets_removed, sample_count, "predicted_celltype_majority_vote", output_dir, naming_token, paste0("_clusters_by_cell_type_majority_vote_", date, "_", i, "_ATAC_overlap.png"))
}




messy_clusters <- c(4, 30) # For multiome F
idxPass <- which(Idents(sc_obj) %in% messy_clusters)
cellsPass <- names(sc_obj$orig.ident[-idxPass])
sc_obj.minus.messy.clusters <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
sc_obj.minus.messy.clusters <- RunPCA(sc_obj.minus.messy.clusters, npcs = 30, approx = T, verbose = T, seed.use = SEED)
sc_obj.minus.messy.clusters <- RunUMAP(sc_obj.minus.messy.clusters, reduction = "pca", dims = 1:30, seed.use = SEED, return.model = T)
sc_obj.minus.messy.clusters <- MajorityVote(sc_obj.minus.messy.clusters)
print_UMAP(sc_obj.minus.messy.clusters, sample_count, "predicted_celltype_majority_vote", output_dir, naming_token, paste0("_clusters_by_cell_type_majority_vote_without_messy_clusters_", date, ".png"))

for (res in seq(0, 3, 0.3)) {
  sc_obj.minus.messy.clusters <- FindClusters(sc_obj.minus.messy.clusters, resolution = res)
}
clustree(sc_obj.minus.messy.clusters, prefix = "integrated_snn_res.")
ggsave(paste0(output_dir, naming_token, "_cluster.trees", date, ".png"), device = "png", width = 8, height = 8, units = "in")




# Code to generate table of which cells are removed at each doublet threshold
removed_doublets_by_cell_type_df <- data.frame("Cell type" = unique(sc_obj$predicted.id))
previous_i <- 1
for (i in seq(0.95, 0.05, by=-0.05)) {
  sc_obj_doublets_removed <- subset(x = sc_obj, subset = scDblFinder.score > i & scDblFinder.score <= previous_i)
  cell_counts <- c()
  for(current_cell_type in unique(sc_obj$predicted.id)) {
    # Grab cells associated with cell type
    idxPass <- which(sc_obj_doublets_removed$predicted.id %in% current_cell_type)
    cellsPass <- names(sc_obj_doublets_removed$orig.ident[idxPass])
    if (length(cellsPass) == 0) {
      cell_counts <- c(cell_counts, 0)
    } else {
      cell_counts <- c(cell_counts, length(cellsPass))
    }
  }
  temp_df <- data.frame(cell_counts)
  names(temp_df)[names(temp_df) == "cell_counts"] <- paste0("Range: ", i, " to ", previous_i)
  removed_doublets_by_cell_type_df <- cbind(removed_doublets_by_cell_type_df, temp_df)
  previous_i <- i
}
summed_cells_removed <- apply(removed_doublets_by_cell_type_df[,2:ncol(removed_doublets_by_cell_type_df)], 1, function(x){sum(x)})
removed_doublets_by_cell_type_df$`Total Cells Removed` <- summed_cells_removed
# Calculate percentage of cells removed per cell type
total_cells <- c()
percentage_of_cells_removed <- c()
i <- 1
for(current_cell_type in unique(sc_obj$predicted.id)) {
  # Grab cells associated with cell type
  idxPass <- which(sc_obj$predicted.id %in% current_cell_type)
  cellsPass <- names(sc_obj$orig.ident[idxPass])
  total_cells <- c(total_cells, length(cellsPass))
  percentage_of_cells_removed <- c(percentage_of_cells_removed, summed_cells_removed[i] / length(cellsPass))
  i <- i + 1
}
removed_doublets_by_cell_type_df$`Total Cells` <- total_cells
removed_doublets_by_cell_type_df$`Percentage of Cells Removed` <- percentage_of_cells_removed
column_sums <- apply(removed_doublets_by_cell_type_df[,2:ncol(removed_doublets_by_cell_type_df)], 2, function(x){sum(x)})
column_sums <- append(column_sums, "Sum", after = 0)
column_sums[length(column_sums)] <- "N/A"
removed_doublets_by_cell_type_df[nrow(removed_doublets_by_cell_type_df) + 1,] = column_sums

write.csv(removed_doublets_by_cell_type_df, file = paste0(output_dir, naming_token, "_removed_doublets_by_cell_type_", date, ".csv"), quote = FALSE, row.names = FALSE)

# Code to generate table of number of doublets removed at each doublet threshold per sample
removed_doublets_by_sample_df <- data.frame("Sample" = unique(sc_obj$sample))
previous_i <- 1
for (i in seq(0.95, 0.05, by=-0.05)) {
  sc_obj_doublets_removed <- subset(x = sc_obj, subset = scDblFinder.score > i & scDblFinder.score <= previous_i)
  sample_counts <- c()
  for(current_sample in unique(sc_obj$sample)) {
    # Grab cells associated with sample
    idxPass <- which(sc_obj_doublets_removed$sample %in% current_sample)
    cellsPass <- names(sc_obj_doublets_removed$orig.ident[idxPass])
    if (length(cellsPass) == 0) {
      sample_counts <- c(sample_counts, 0)
    } else {
      sample_counts <- c(sample_counts, length(cellsPass))
    }
  }
  temp_df <- data.frame(sample_counts)
  names(temp_df)[names(temp_df) == "sample_counts"] <- paste0("Range: ", i, " to ", previous_i)
  removed_doublets_by_sample_df <- cbind(removed_doublets_by_sample_df, temp_df)
  previous_i <- i
}
summed_cells_removed <- apply(removed_doublets_by_sample_df[,2:ncol(removed_doublets_by_sample_df)], 1, function(x){sum(x)})
removed_doublets_by_sample_df$`Total Cells Removed` <- summed_cells_removed
# Calculate percentage of cells removed per cell type
total_cells <- c()
percentage_of_cells_removed <- c()
i <- 1
for(current_sample in unique(sc_obj$sample)) {
  # Grab cells associated with cell type
  idxPass <- which(sc_obj$sample %in% current_sample)
  cellsPass <- names(sc_obj$orig.ident[idxPass])
  total_cells <- c(total_cells, length(cellsPass))
  percentage_of_cells_removed <- c(percentage_of_cells_removed, summed_cells_removed[i] / length(cellsPass))
  i <- i + 1
}
removed_doublets_by_sample_df$`Total Cells` <- total_cells
removed_doublets_by_sample_df$`Percentage of Cells for Cell Type Removed` <- percentage_of_cells_removed
column_sums <- apply(removed_doublets_by_sample_df[,2:ncol(removed_doublets_by_sample_df)], 2, function(x){sum(x)})
column_sums <- append(column_sums, "Sum", after = 0)
column_sums[length(column_sums)] <- "N/A"
removed_doublets_by_sample_df[nrow(removed_doublets_by_sample_df) + 1,] = column_sums

write.csv(removed_doublets_by_sample_df, file = paste0(output_dir, naming_token, "_removed_doublets_by_sample_", date, ".csv"), quote = FALSE, row.names = FALSE)







# Feature plot on doublet score
FeaturePlot(sc_obj, features = "scDblFinder.score")
ggsave(paste0(output_dir, naming_token, "_doublet_score_feature_plot_", date, ".png"), device = "png", dpi = 300, width = 20, height = 20, units = "in")

# Ridge plot on doublet score (per sample)
ridge_plot_for_doublets <- RidgePlot(sc_obj, features = 'scDblFinder.score', group.by = "sample")
ridge_plot_for_doublets <- ridge_plot_for_doublets + scale_x_continuous(breaks = scales::pretty_breaks(n = 20))
ggsave(paste0(output_dir, naming_token, "_doublet_score_ridge_plot_", date, ".png"), device = "png", dpi = 300, width = 20, height = 20, units = "in")


# Code to parse marker tables
marker_dir <- "C:/Users/willi/Desktop/multiome junk/02-24/combined_cell_types/2.4/cluster_markers/SCT"
output_dir <- "C:/Users/willi/Desktop/multiome junk/02-24/combined_cell_types/2.4/cluster_markers/SCT/UPDATED/"
marker_files <- list.files(marker_dir, pattern = "*.txt$", full.names = TRUE)
for(marker_file in marker_files) {
  marker_file_content <- read.table(marker_file)
  marker_file_content <- marker_file_content[marker_file_content$avg_log2FC > 0,]
  marker_file_content <- marker_file_content[marker_file_content$p_val_adj < 0.05,]
  marker_file_content <- marker_file_content[order(marker_file_content$pct.1, decreasing = TRUE),]
  write.table(marker_file_content, paste0(output_dir, basename(marker_file)), quote = FALSE, sep = "\t")
}


# Code to parse marker table in one file
marker_file <- "C:/Users/wat2/Desktop/02-28/7_high_vs_low_viral_load_D28_V3_cluster_markers_scaled.txt"
output_dir <- "C:/Users/wat2/Desktop/02-28/UPDATED/"
marker_file_content <- read.table(marker_file)
for(cluster_id in unique(marker_file_content$cluster)) {
  current_cluster_content <- marker_file_content[marker_file_content$cluster == cluster_id,]
  current_cluster_content <- current_cluster_content[current_cluster_content$avg_log2FC > 0,]
  current_cluster_content <- current_cluster_content[current_cluster_content$p_val_adj < 0.05,]
  current_cluster_content <- current_cluster_content[order(current_cluster_content$pct.1, decreasing = TRUE),]
  write.table(current_cluster_content, paste0(output_dir, cluster_id, "_", basename(marker_file)), quote = FALSE, sep = "\t")
}

# Code to parse DEG markers into positive and negative fold change and change Bonferroni to FDR
marker_dir <- "C:/Users/wat2/Desktop/02-28/DEGs/Unfiltered"
output_dir <- "C:/Users/wat2/Desktop/02-28/DEGs/Unfiltered/UPDATED/"
marker_files <- list.files(marker_dir, pattern = "*.csv$", full.names = TRUE)
for(marker_file in marker_files) {
  marker_file_content <- read.table(marker_file, sep = ",", header = TRUE)
  marker_file_content$p_val_adj <- p.adjust(marker_file_content$p_val, method='fdr')
  positive_fc_markers <- marker_file_content[marker_file_content$avg_log2FC > 0,]
  negative_fc_markers <- marker_file_content[marker_file_content$avg_log2FC < 0,]
  write.table(positive_fc_markers, paste0(output_dir, "UNFILTERED_POS_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  write.table(negative_fc_markers, paste0(output_dir, "UNFILTERED_NEG_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  # Should I adjust p value based on only positive or negative values? Probably not
  positive_fc_markers <- positive_fc_markers[positive_fc_markers$p_val_adj < 0.05,]
  negative_fc_markers <- negative_fc_markers[negative_fc_markers$p_val_adj < 0.05,]
  write.table(positive_fc_markers, paste0(output_dir, "POS_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  write.table(positive_fc_markers$X, paste0(output_dir, "GENES_POS_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
  write.table(negative_fc_markers, paste0(output_dir, "NEG_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  write.table(negative_fc_markers$X, paste0(output_dir, "GENES_NEG_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
  # 0.1, 0.585, 1, 2
  positive_fc_markers_0.1 <- positive_fc_markers[positive_fc_markers$avg_log2FC > 0.1,]
  positive_fc_markers_0.585 <- positive_fc_markers[positive_fc_markers$avg_log2FC > 0.585,]
  positive_fc_markers_1 <- positive_fc_markers[positive_fc_markers$avg_log2FC > 1,]
  positive_fc_markers_2 <- positive_fc_markers[positive_fc_markers$avg_log2FC > 2,]
  negative_fc_markers_0.1 <- negative_fc_markers[negative_fc_markers$avg_log2FC < -0.1,]
  negative_fc_markers_0.585 <- negative_fc_markers[negative_fc_markers$avg_log2FC < -0.585,]
  negative_fc_markers_1 <- negative_fc_markers[negative_fc_markers$avg_log2FC < -1,] 
  negative_fc_markers_2 <- negative_fc_markers[negative_fc_markers$avg_log2FC < -2,]
  write.table(positive_fc_markers_0.1, paste0(output_dir, "0.1_POS_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  
  write.table(positive_fc_markers_0.1$X, paste0(output_dir, "GENES_0.1_POS_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  
  write.table(positive_fc_markers_0.585, paste0(output_dir, "0.585_POS_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  
  write.table(positive_fc_markers_0.585$X, paste0(output_dir, "GENES_0.585_POS_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  
  write.table(positive_fc_markers_1, paste0(output_dir, "1_POS_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  
  write.table(positive_fc_markers_1$X, paste0(output_dir, "GENES_1_POS_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  
  write.table(positive_fc_markers_2, paste0(output_dir, "2_POS_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  
  write.table(positive_fc_markers_2$X, paste0(output_dir, "GENES_2_POS_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  
  write.table(negative_fc_markers_0.1, paste0(output_dir, "0.1_NEG_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  
  write.table(negative_fc_markers_0.1$X, paste0(output_dir, "GENES_0.1_NEG_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  
  write.table(negative_fc_markers_0.585, paste0(output_dir, "0.585_NEG_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  
  write.table(negative_fc_markers_0.585$X, paste0(output_dir, "GENES_0.585_NEG_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  
  write.table(negative_fc_markers_1, paste0(output_dir, "1_NEG_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  
  write.table(negative_fc_markers_1$X, paste0(output_dir, "GENES_1_NEG_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
  
  write.table(negative_fc_markers_2, paste0(output_dir, "2_NEG_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)

  write.table(negative_fc_markers_2$X, paste0(output_dir, "GENES_2_NEG_", basename(marker_file)), quote = FALSE, sep = ",", row.names = FALSE)
}


