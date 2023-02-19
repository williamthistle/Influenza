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
  print_UMAP(sc_obj, sample_count, "predicted_celltype_majority_vote", output_dir, naming_token, paste0("_clusters_by_cell_type_majority_vote_", date, SEED, ".png"))
  print_UMAP(sc_obj, sample_count, "seurat_clusters", output_dir, naming_token, paste0("_clusters_by_cluster_num_", date, SEED, ".png"))
  print_UMAP(sc_obj, sample_count, "predicted.id", output_dir, naming_token, paste0("_clusters_by_cell_type_", date, SEED, ".png"))
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
