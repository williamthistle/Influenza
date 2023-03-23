run_differential_expression_cluster <- function(sc_obj, marker_dir) {
  print(paste0("Performing differential expression for each cluster"))
  cluster_ids <- unique(sc_obj$seurat_clusters)
  for(cluster_id in cluster_ids) {
    print(cluster_id)
    cluster.markers <- FindMarkers(sc_obj, ident.1 = cluster_id, assay = "SCT")
    write.table(cluster.markers, paste0(marker_dir, "7_", cluster_id, ".txt"), quote = FALSE, sep = "\t")
  }
}

run_differential_expression_group <- function(sc_obj, analysis_dir, group) {
  print(paste0("Performing differential expression for group ", group, " for each cell type"))
  
  # Make reading data parallel
  if (Sys.getenv("SLURM_NTASKS_PER_NODE") == "") {
    n.cores <- as.numeric(detectCores())
  } else {
    n.cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
  }
  
  n.cores <- 8
  
  message(paste0("Number of cores: ", n.cores))
  
  registerDoMC(n.cores)
  message("Begin parallelizing...")
  
  all_cell_types <- union(unique(sc_obj$predicted_celltype_majority_vote), unique(sc_obj$magical_cell_types))
  
  compiled_output <- foreach(
    i = 1:length(all_cell_types),
    .packages = c("Seurat", "base")
  ) %dopar% {
    cell_type <- all_cell_types[[i]]
    print(cell_type)
    if(cell_type %in% unique(sc_obj$predicted_celltype_majority_vote)) {
      idxPass <- which(sc_obj$predicted_celltype_majority_vote %in% cell_type)
    } else {
      idxPass <- which(sc_obj$magical_cell_types %in% cell_type)
      print("This cell type is for MAGICAL processing")
    }
    print(length(idxPass))
    cellsPass <- names(sc_obj$orig.ident[idxPass])
    cells_subset <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
    # TODO: Make this print table for relevant group
    #print(table(cells_subset$viral_load))
    DefaultAssay(cells_subset) <- "SCT"
    Idents(cells_subset) <- group
    if(group == "viral_load") {
      first_group <- "HIGH"
      second_group <- "LOW"
    } else if(group == "day") {
      first_group <- "D28"
      second_group <- "D_MINUS_1"
    } else if(group == "sex") {
      first_group <- "MALE"
      second_group <- "FEMALE"
    }
    diff_markers <- FindMarkers(cells_subset, ident.1 = first_group, ident.2 = second_group, assay = "SCT", recorrect_umi = FALSE, logfc.threshold = 0, min.pct = 0)
    cell_type <- sub(" ", "_", cell_type)
    write.csv(diff_markers, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type, "-", group, ".csv"), quote = FALSE)
    return(i)
  }
  print(paste0("Done performing differential expression for group ", group, " for each cell type"))
}