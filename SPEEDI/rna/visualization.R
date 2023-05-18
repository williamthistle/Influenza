# Method to print a well-organized UMAP plot for our snRNA-seq data
print_UMAP <- function(sc_obj, sample_count, group_by_category, plot_dir, file_name) {
  cell_count <- length(sc_obj$cell_name)
  current_title <- paste0("snRNA-seq Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
  DimPlot(sc_obj, reduction = "umap", group.by = group_by_category, label = TRUE,
          label.size = 3, repel = TRUE, raster = FALSE) +
    labs(title = current_title) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(plot_dir, file_name), device = "png", dpi = 300)
}

print_UMAP_stage_1 <- function(sc_obj, sample_count, plot_dir, date) {
  print_UMAP(sc_obj, sample_count, "predicted_celltype_majority_vote", plot_dir, paste0("pre.clusters_by_cell_type_majority_vote_", date, ".png"))
  print_UMAP(sc_obj, sample_count, "seurat_clusters", plot_dir, paste0("pre.clusters_by_cluster_num_", date, ".png"))
  print_UMAP(sc_obj, sample_count, "predicted.id", plot_dir, paste0("pre.clusters_by_cell_type_", date, ".png"))
  print_UMAP(sc_obj, sample_count, "viral_load", plot_dir, paste0("pre.clusters_by_viral_load_", date, ".png"))
  print_UMAP(sc_obj, sample_count, "sample", plot_dir, paste0("pre.clusters_by_sample_", date, ".png"))
}

print_UMAP_tagged <- function(sc_obj, tagged_sample_list, sample_count, plot_dir, date) {
  # Tag certain samples to see whether they're very different from other samples
  sample_vec <- sc_obj$sample
  sample_vec[!(sample_vec %in% tagged_sample_list)] <- "UNTAGGED"
  sample_vec[sample_vec %in% tagged_sample_list] <- "TAGGED"
  sc_obj$tagged <- sample_vec
  print_UMAP(sc_obj, sample_count, "tagged", plot_dir, paste0("pre.clusters_by_tag_", date, ".png"))
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