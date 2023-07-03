# Print UMAP (RNA)
print_UMAP_RNA <- function(sc_obj, file_name, group_by_category = NULL, output_dir = getwd(), log_flag = FALSE) {
  # Normalize paths (in case user provides relative paths)
  output_dir <- normalize_dir_path(output_dir)
  sample_count <- length(unique(sc_obj$sample))
  cell_count <- length(sc_obj$cell_name)
  current_title <- paste0("RNA Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
  if(!is.null(group_by_category)) {
    p <- Seurat::DimPlot(sc_obj, reduction = "umap", group.by = group_by_category, label = TRUE,
                         label.size = 3, repel = TRUE, raster = FALSE) +
      ggplot2::labs(title = current_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  } else {
    p <- Seurat::DimPlot(sc_obj, reduction = "umap", label = TRUE,
                         label.size = 3, repel = TRUE, raster = FALSE) +
      ggplot2::labs(title = current_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  ggplot2::ggsave(paste0(output_dir, file_name), plot = p, device = "png", dpi = 300)
  return(TRUE)
}

# Initial combination of cell types
combine_cell_types_initial <- function(sc_obj, resolution = 1.5) {
  sc_obj$old.predicted.id <- sc_obj$predicted.id
  Cell_type_combined <- sc_obj$predicted.id
  idx <- grep("CD4 T", Cell_type_combined)
  Cell_type_combined[idx] <- "CD4 Memory"
  idx <- grep("CD8 T", Cell_type_combined)
  Cell_type_combined[idx] <- "CD8 Memory"
  idx <- grep("cDC", Cell_type_combined)
  Cell_type_combined[idx] <- "cDC"
  idx <- grep("Proliferating", Cell_type_combined)
  Cell_type_combined[idx] <- "Proliferating"
  sc_obj$predicted.id <- Cell_type_combined
  sc_obj <- MajorityVote(sc_obj, resolution)
  return(sc_obj)
}

# Combine more cell types (for MAGICAL)
combine_cell_types_magical <- function(sc_obj) {
  Cell_type_combined <- sc_obj$predicted_celltype_majority_vote
  levels(Cell_type_combined) <- c(levels(Cell_type_combined), "T Naive", "B", "NK_MAGICAL")
  idx <- grep("CD4 Naive", Cell_type_combined)
  Cell_type_combined[idx] <- "T Naive"
  idx <- grep("CD8 Naive", Cell_type_combined)
  Cell_type_combined[idx] <- "T Naive"
  idx <- grep("Treg", Cell_type_combined)
  Cell_type_combined[idx] <- "T Naive"
  idx <- grep("NK", Cell_type_combined)
  Cell_type_combined[idx] <- "NK_MAGICAL"
  idx <- grep("B", Cell_type_combined)
  Cell_type_combined[idx] <- "B"
  sc_obj$magical_cell_types <- Cell_type_combined
  return(sc_obj)
}

# Override cluster label from majority voting
override_cluster_label <- function(sc_obj, cluster_nums, new_cluster_label) {
  for(cluster_id in cluster_nums) {
    print(cluster_id)
    if(!(new_cluster_label %in% levels(sc_obj$predicted_celltype_majority_vote))) {
      levels(sc_obj$predicted_celltype_majority_vote) <- c(levels(sc_obj$predicted_celltype_majority_vote), new_cluster_label)
    }
    majority_vote <- sc_obj$predicted_celltype_majority_vote
    idxPass <- which(Idents(sc_obj) %in% c(cluster_id))
    majority_vote[idxPass] <- new_cluster_label
    sc_obj$predicted_celltype_majority_vote <- majority_vote
  }
  return(sc_obj)
}

# Remove specific samples from Seurat object
remove_specific_samples_from_sc_obj <- function(sc_obj, samples) {
  idxPass <- which(sc_obj$sample %in% samples)
  cellsPass <- names(sc_obj$orig.ident[-idxPass])
  sc_obj <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
  return(sc_obj)
}

remove_cells_based_on_umap <- function(sc_obj, first_x, second_x, first_y, second_y) {
  orig.umap.coords <- as.data.frame(sc_obj[["umap"]]@cell.embeddings)
  orig.umap.coords$cells <- rownames(orig.umap.coords)
  deleted.cells.umap.coords <- orig.umap.coords[orig.umap.coords$"UMAP_1" > first_x & orig.umap.coords$"UMAP_1" < second_x,]
  deleted.cells.umap.coords <- deleted.cells.umap.coords[deleted.cells.umap.coords$"UMAP_2" > first_y & deleted.cells.umap.coords$"UMAP_2" < second_y,]
  final.umap.coords <- orig.umap.coords[!(orig.umap.coords$cells %in% deleted.cells.umap.coords$cells),]
  cellsPass <- rownames(final.umap.coords)
  sc_obj <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
  return(sc_obj)
}

run_differential_expression_cluster <- function(sc_obj, marker_dir) {
  print(paste0("Performing differential expression for each cluster"))
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
  cluster_ids <- unique(sc_obj$seurat_clusters)
  
  compiled_output <- foreach(
    i = 1:length(cluster_ids),
    .packages = c("Seurat", "base")
  ) %dopar% {
    cluster_id <- cluster_ids[[i]]
    cluster.markers <- FindMarkers(sc_obj, ident.1 = cluster_id, assay = "RNA")
    write.table(cluster.markers, paste0(marker_dir, "7_", cluster_id, ".txt"), quote = FALSE, sep = "\t")
    return(i)
  }
  message("All done!")
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
    DefaultAssay(cells_subset) <- "RNA"
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
    #diff_markers <- FindMarkers(cells_subset, ident.1 = first_group, ident.2 = second_group, logfc.threshold = 0, min.pct = 0)
    diff_markers = FindMarkers(cells_subset, test.use="LR", latent.vars = 'subject', ident.1 = first_group, ident.2 = second_group, logfc.threshold = 0, min.pct = 0)
    cell_type <- sub(" ", "_", cell_type)
    write.csv(diff_markers, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type, "-", group, ".csv"), quote = FALSE)
    return(i)
  }
  print(paste0("Done performing differential expression for group ", group, " for each cell type"))
}