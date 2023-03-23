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