# Method to print a well-organized UMAP plot for our snRNA-seq data
print_UMAP_RNA_for_paper <- function(sc_obj, file_name, group_by_category = NULL, output_dir = getwd(), log_flag = FALSE) {
  # Normalize paths (in case user provides relative paths)
  output_dir <- normalize_dir_path(output_dir)
  sample_count <- length(unique(sc_obj$sample))
  cell_count <- length(sc_obj$cell_name)
  current_title <- paste0("RNA Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
  if(!is.null(group_by_category)) {
    p <- Seurat::DimPlot(sc_obj, reduction = "umap", group.by = group_by_category, repel = TRUE, raster = FALSE) + 
      labs(x ="UMAP 1", y = "UMAP 2", title = NULL) + theme_classic(base_size = 24)
    
  } else {
    p <- Seurat::DimPlot(sc_obj, reduction = "umap", repel = TRUE, raster = FALSE) +
      labs(x ="UMAP 1", y = "UMAP 2", title = NULL) + theme_classic(base_size = 24)
  }
  ggplot2::ggsave(paste0(output_dir, file_name), plot = p, device = "png", dpi = 300, width=8.5, height=5.5, units = "in")
  return(TRUE)
}

# Method to print a well-organized UMAP plot for our snRNA-seq data
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

# Method to print a well-organized UMAP plot for our snRNA-seq data (with extra axis ticks)
print_UMAP_RNA_detailed <- function(sc_obj, file_name, group_by_category = NULL, output_dir = getwd(), log_flag = FALSE) {
  # Normalize paths (in case user provides relative paths)
  output_dir <- normalize_dir_path(output_dir)
  sample_count <- length(unique(sc_obj$sample))
  cell_count <- length(sc_obj$cell_name)
  current_title <- paste0("RNA Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
  if(!is.null(group_by_category)) {
    p <- Seurat::DimPlot(sc_obj, reduction = "umap", group.by = group_by_category, label = TRUE,
                         label.size = 3, repel = TRUE, raster = FALSE) +
      ggplot2::labs(title = current_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + scale_x_continuous(n.breaks=20) +
      scale_y_continuous(n.breaks=20)
  } else {
    p <- Seurat::DimPlot(sc_obj, reduction = "umap", label = TRUE,
                         label.size = 3, repel = TRUE, raster = FALSE) +
      ggplot2::labs(title = current_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + scale_x_continuous(n.breaks=20) +
      scale_y_continuous(n.breaks=20)
  }
  ggplot2::ggsave(paste0(output_dir, file_name), plot = p, device = "png", dpi = 300)
  return(TRUE)
}

print_UMAP_stage_1 <- function(sc_obj, sample_count, plot_dir, date) {
  print_UMAP(sc_obj, sample_count, "predicted_celltype_majority_vote", plot_dir, paste0("pre.clusters_by_cell_type_majority_vote_", date, ".png"))
  print_UMAP(sc_obj, sample_count, "seurat_clusters", plot_dir, paste0("pre.clusters_by_cluster_num_", date, ".png"))
  print_UMAP(sc_obj, sample_count, "predicted.id", plot_dir, paste0("pre.clusters_by_cell_type_", date, ".png"))
  print_UMAP(sc_obj, sample_count, "viral_load", plot_dir, paste0("pre.clusters_by_viral_load_", date, ".png"))
  print_UMAP(sc_obj, sample_count, "sample", plot_dir, paste0("pre.clusters_by_sample_", date, ".png"))
}

print_final_UMAPs <- function(sc_obj, output_dir, token = "Default") {
  print_UMAP_RNA(sc_obj, file_name = paste0(token, "_RNA_UMAP_by_Majority_Vote_Cell_Type.png"),
                 group_by_category = "predicted_celltype_majority_vote", output_dir = RNA_output_dir)
  print_UMAP_RNA(sc_obj, file_name = paste0(token, "_RNA_UMAP_by_Cluster.png"),
                 group_by_category = "seurat_clusters", output_dir = RNA_output_dir)
  print_UMAP_RNA(sc_obj, file_name = paste0(token, "_RNA_UMAP_by_Raw_Predicted_Cell_Type.png"),
                 group_by_category = "predicted.id", output_dir = RNA_output_dir)
  print_UMAP_RNA(sc_obj, file_name = paste0(token, "_RNA_UMAP_by_Sample.png"),
                 group_by_category = "sample", output_dir = RNA_output_dir)
  print_UMAP_RNA(sc_obj, file_name = paste0(token, "_RNA_UMAP_by_Day.png"),
                 group_by_category = "time_point", output_dir = RNA_output_dir)
  print_UMAP_RNA(sc_obj, file_name = paste0(token, "_RNA_UMAP_by_Sex.png"),
                 group_by_category = "sex", output_dir = RNA_output_dir)
}

print_pheatmap_predicted_id <- function(sc_obj, plot_dir, date, all_viral_load_samples, high_viral_load_samples, low_viral_load_samples) {
  # Let's also use Vincy's code to create a nice heatmap of cell type proportions
  sc_obj$sample <- factor(sc_obj$sample, levels = all_viral_load_samples)
  voting_cell_type_proportion <- as.matrix(table(sc_obj$sample, sc_obj$predicted.id))
  voting_cell_type_proportion <- apply(voting_cell_type_proportion, 1, function(x){x/sum(x)})
  my_sample_col <- data.frame(sample = rep(c("HIGH", "LOW"), c(length(high_viral_load_samples),length(low_viral_load_samples)))) #The study group by their order
  row.names(my_sample_col) <- colnames(voting_cell_type_proportion)
  output.plot <- pheatmap(voting_cell_type_proportion, annotation_col = my_sample_col, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.3f")
  ggsave(paste0(plot_dir, "voting_cell_type_proportions_heatmap_predicted.id_", date, ".png"), plot = output.plot, device = "png", width = 8, height = 8, units = "in")
}

print_pheatmap_majority_vote <- function(sc_obj, plot_dir, date, all_viral_load_samples, high_viral_load_samples, low_viral_load_samples) {
  # Let's also use Vincy's code to create a nice heatmap of cell type proportions
  sc_obj$sample <- factor(sc_obj$sample, levels = all_viral_load_samples)
  voting_cell_type_proportion <- as.matrix(table(sc_obj$sample, sc_obj$predicted_celltype_majority_vote))
  voting_cell_type_proportion <- apply(voting_cell_type_proportion, 1, function(x){x/sum(x)})
  my_sample_col <- data.frame(sample = rep(c("HIGH", "LOW"), c(length(high_viral_load_samples),length(low_viral_load_samples)))) #The study group by their order
  row.names(my_sample_col) <- colnames(voting_cell_type_proportion)
  output.plot <- pheatmap(voting_cell_type_proportion, annotation_col = my_sample_col, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.3f")
  ggsave(paste0(plot_dir, "voting_cell_type_proportions_heatmap_majority_vote_", date, ".png"), plot = output.plot, device = "png", width = 8, height = 8, units = "in")
}