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
