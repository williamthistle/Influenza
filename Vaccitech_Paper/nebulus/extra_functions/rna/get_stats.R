# Print clustree plot to figure out best resolution for clustering
print_clustree_plot <- function(sc_obj, plot_dir, date) {
  for (res in seq(0, 6, 0.3)) {
    sc_obj <- FindClusters(sc_obj, resolution = res)
  }
  clustree(sc_obj, prefix = "integrated_snn_res.")
  ggsave(paste0(plot_dir, "cluster.trees_", date, ".png"), device = "png", width = 20, height = 20, units = "in")
}

# Check how many cells from each condition are present for each cell type
print_celltype_counts <- function(sc_obj) {
  for (cell_type in unique(sc_obj$predicted_celltype_majority_vote)) {
    print(cell_type)
    idxPass <- which(sc_obj$predicted_celltype_majority_vote %in% cell_type)
    cellsPass <- names(sc_obj$orig.ident[idxPass])
    sample_subset <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
    print(length(cellsPass))
    print(table(sample_subset$viral_load))
    print(table(sample_subset$day))
    print(table(sample_subset$sex))
  }
}

# Calculate cell type proportions and counts (not really being used currently)
calculate_props_and_counts <- function(sc_obj, condition, sample_names) {
  cell_type_proportions_df <- data.frame("Condition" = condition, "Sample_name" = sample_names)
  total_cell_counts_df <- data.frame("Condition" = condition, "Sample_name" = sample_names)
  total_cell_counts <- c()
  for (sample_id in sample_names) {
    idxPass <- which(sc_obj$sample %in% sample_id)
    cellsPass <- names(sc_obj$orig.ident[idxPass])
    sample_subset <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
    total_cell_counts <- c(total_cell_counts, ncol(sample_subset))
  }
  total_cell_counts_df <- cbind(total_cell_counts_df, total_cell_counts)
  for (cell_type in unique(sc_obj$predicted.id)) {
    cell_type_proportions <- vector()
    cell_counts <- c()
    print(cell_type)
    # Grab cells associated with cell type
    idxPass <- which(sc_obj$predicted.id %in% cell_type)
    print(length(idxPass))
    cellsPass <- names(sc_obj$orig.ident[idxPass])
    cells_subset <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
    for (sample_id in sample_names) {
      # Subset further based on cells associated with sample ID
      idxPass <- which(cells_subset$sample %in% sample_id)
      cellsPass <- names(cells_subset$orig.ident[idxPass])
      if (length(cellsPass) == 0) {
        cell_counts <- c(cell_counts, 0)
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
  return(list(cell_type_proportions_df, total_cell_counts_df))
}

# Capture info about clusters
capture_cluster_info <- function(sc_obj, find_doublet_info = FALSE) {
  cluster_ids <- vector()
  cluster_distributions <- list()
  cluster_predictions <- vector()
  cluster_mean_S_score <- vector()
  cluster_mean_G2M_score <- vector()
  cluster_mean_CC_difference <- vector()
  num_cells <- c()
  num_cells_high <- c()
  num_cells_low <- c()
  cluster_mean_mito <- c()
  cluster_mean_mito_high <- c()
  cluster_mean_mito_low <- c()
  cluster_mean_nFeature <- c()
  cluster_mean_nFeature_high <- c()
  cluster_mean_nFeature_low <- c()
  cluster_mean_nCount <- c()
  cluster_mean_nCount_high <- c()
  cluster_mean_nCount_low <- c()
  cluster_mean_rp <- c()
  cluster_mean_rp_high <- c()
  cluster_mean_rp_low <- c()
  cluster_mean_doublet <- c()
  cluster_mean_doublet_high <- c()
  cluster_mean_doublet_low <- c()
  idx <- 1
  for (cluster in levels(sc_obj)) {
    cluster_ids <- append(cluster_ids, cluster)
    idxPass <- which(Idents(sc_obj) %in% cluster)
    cellsPass <- names(sc_obj$orig.ident[idxPass])
    num_cells <- c(num_cells, length(cellsPass))
    filtered_cluster <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
    cluster_mean_mito <- c(cluster_mean_mito, mean(filtered_cluster$percent.mt))
    cluster_mean_nFeature <- c(cluster_mean_nFeature, mean(filtered_cluster$nFeature_RNA))
    cluster_mean_nCount <- c(cluster_mean_nCount, mean(filtered_cluster$nCount_RNA))
    cluster_mean_rp <- c(cluster_mean_rp, mean(filtered_cluster$percent.rp))
    if(find_doublet_info) {
      cluster_mean_doublet <- c(cluster_mean_doublet, mean(filtered_cluster$scDblFinder.score))
    }
    cluster_prediction <- sort(table(filtered_cluster$predicted_celltype_majority_vote), decreasing = TRUE)[1]
    cluster_predictions <- append(cluster_predictions, cluster_prediction)
    cluster_distributions[[idx]] <- table(filtered_cluster$predicted.id)
    cluster_mean_S_score <- append(cluster_mean_S_score, mean(filtered_cluster$S.Score))
    cluster_mean_G2M_score <- append(cluster_mean_G2M_score, mean(filtered_cluster$G2M.Score))
    cluster_mean_CC_difference <- append(cluster_mean_CC_difference, mean(filtered_cluster$CC.Difference))
    # High
    idxPass <- which(filtered_cluster$viral_load_category %in% "high")
    cellsPass <- names(filtered_cluster$orig.ident[idxPass])
    num_cells_high <- c(num_cells_high, length(cellsPass))
    if(length(cellsPass) > 100) {
      filtered_cluster_high <- subset(x = filtered_cluster, subset = cell_name %in% cellsPass)
      cluster_mean_mito_high <- c(cluster_mean_mito_high, mean(filtered_cluster_high$percent.mt))
      cluster_mean_nFeature_high <- c(cluster_mean_nFeature_high, mean(filtered_cluster_high$nFeature_RNA))
      cluster_mean_nCount_high <- c(cluster_mean_nCount_high, mean(filtered_cluster_high$nCount_RNA))
      cluster_mean_rp_high <- c(cluster_mean_rp_high, mean(filtered_cluster_high$percent.rp))
      if(find_doublet_info) {
        cluster_mean_doublet_high <- c(cluster_mean_doublet_high, mean(filtered_cluster_high$scDblFinder.score))
      }
    } else {
      cluster_mean_mito_high <- c(cluster_mean_mito_high, NA)
      cluster_mean_nFeature_high <- c(cluster_mean_nFeature_high, NA)
      cluster_mean_nCount_high <- c(cluster_mean_nCount_high, NA)
      cluster_mean_rp_high <- c(cluster_mean_rp_high, NA)
      if(find_doublet_info) {
        cluster_mean_doublet_high <- c(cluster_mean_doublet_high, NA)
      }
    }
    # Low
    idxPass <- which(filtered_cluster$viral_load_category %in% "low")
    cellsPass <- names(filtered_cluster$orig.ident[idxPass])
    num_cells_low <- c(num_cells_low, length(cellsPass))
    if(length(cellsPass) > 100) {
      filtered_cluster_low <- subset(x = filtered_cluster, subset = cell_name %in% cellsPass)
      cluster_mean_mito_low <- c(cluster_mean_mito_low, mean(filtered_cluster_low$percent.mt))
      cluster_mean_nFeature_low <- c(cluster_mean_nFeature_low, mean(filtered_cluster_low$nFeature_RNA))
      cluster_mean_nCount_low <- c(cluster_mean_nCount_low, mean(filtered_cluster_low$nCount_RNA))
      cluster_mean_rp_low <- c(cluster_mean_rp_low, mean(filtered_cluster_low$percent.rp))
      if(find_doublet_info) {
        cluster_mean_doublet_low <- c(cluster_mean_doublet_low, mean(filtered_cluster_low$scDblFinder.score))
      }
    } else {
      cluster_mean_mito_low <- c(cluster_mean_mito_low, NA)
      cluster_mean_nFeature_low <- c(cluster_mean_nFeature_low, NA)
      cluster_mean_nCount_low <- c(cluster_mean_nCount_low, NA)
      cluster_mean_rp_low <- c(cluster_mean_rp_low, NA)
      if(find_doublet_info) {
        cluster_mean_doublet_low <- c(cluster_mean_doublet_low, NA)
      }
    }
    idx <- idx + 1
  }
  names(cluster_predictions) <- paste(levels(sc_obj), "-", names(cluster_predictions))
  cell_cycle_df <- data.frame("Cluster" = cluster_ids, "S" = cluster_mean_S_score, "G2M" = cluster_mean_G2M_score, "CC Diff" = cluster_mean_CC_difference)
  cluster_QC_stats_sorted_by_mt <- data.frame("cluster" = cluster_ids, "num_cells" = num_cells, "num_cells_high_viral" = num_cells_high, "num_cells_low_viral" = num_cells_low,
                                              "mean_mito" = cluster_mean_mito, "mean_mito_high_viral" = cluster_mean_mito_high, "mean_mito_low_viral" = cluster_mean_mito_low,
                                              "mean_mito_viral_diff" = abs(cluster_mean_mito_high - cluster_mean_mito_low))
  rownames(cluster_QC_stats_sorted_by_mt) <- cluster_QC_stats_sorted_by_mt$cluster
  cluster_QC_stats_sorted_by_mt <- cluster_QC_stats_sorted_by_mt[ , !(names(cluster_QC_stats_sorted_by_mt) %in% c("cluster"))]
  cluster_QC_stats_sorted_by_mt <- cluster_QC_stats_sorted_by_mt[order(cluster_QC_stats_sorted_by_mt$mean_mito, decreasing = TRUE),]
  cluster_QC_stats_sorted_by_nFeature <- data.frame("cluster" = cluster_ids, "num_cells" = num_cells, "num_cells_high_viral" = num_cells_high, "num_cells_low_viral" = num_cells_low,
                                                    "mean_nFeature" = cluster_mean_nFeature, "mean_nFeature_high_viral" = cluster_mean_nFeature_high,
                                                    "mean_nFeature_low_viral" = cluster_mean_nFeature_low, "mean_nFeature_diff" = abs(cluster_mean_nFeature_high - cluster_mean_nFeature_low))
  rownames(cluster_QC_stats_sorted_by_nFeature) <- cluster_QC_stats_sorted_by_nFeature$cluster
  cluster_QC_stats_sorted_by_nFeature <- cluster_QC_stats_sorted_by_nFeature[ , !(names(cluster_QC_stats_sorted_by_nFeature) %in% c("cluster"))]
  cluster_QC_stats_sorted_by_nFeature <- cluster_QC_stats_sorted_by_nFeature[order(cluster_QC_stats_sorted_by_nFeature$mean_nFeature, decreasing = TRUE),]
  cluster_QC_stats_sorted_by_nCount <- data.frame("cluster" = cluster_ids, "num_cells" = num_cells, "num_cells_high_viral" = num_cells_high, "num_cells_low_viral" = num_cells_low,
                                                  "mean_nCount" = cluster_mean_nCount, "mean_nCount_high_viral" = cluster_mean_nCount_high,
                                                  "mean_nCount_low_viral" = cluster_mean_nCount_low, "mean_nCount_diff" = abs(cluster_mean_nCount_high - cluster_mean_nCount_low))
  rownames(cluster_QC_stats_sorted_by_nCount) <- cluster_QC_stats_sorted_by_nCount$cluster
  cluster_QC_stats_sorted_by_nCount <- cluster_QC_stats_sorted_by_nCount[ , !(names(cluster_QC_stats_sorted_by_nCount) %in% c("cluster"))]
  cluster_QC_stats_sorted_by_nCount <- cluster_QC_stats_sorted_by_nCount[order(cluster_QC_stats_sorted_by_nCount$mean_nCount, decreasing = TRUE),]
  
  cluster_QC_stats_sorted_by_rp <- data.frame("cluster" = cluster_ids, "num_cells" = num_cells, "num_cells_high_viral" = num_cells_high, "num_cells_low_viral" = num_cells_low,
                                              "mean_rp" = cluster_mean_rp, "mean_rp_high_viral" = cluster_mean_rp_high,
                                              "mean_rp_low_viral" = cluster_mean_rp_low, "mean_rp_diff" = abs(cluster_mean_rp_high - cluster_mean_rp_low))
  rownames(cluster_QC_stats_sorted_by_rp) <- cluster_QC_stats_sorted_by_rp$cluster
  cluster_QC_stats_sorted_by_rp <- cluster_QC_stats_sorted_by_rp[ , !(names(cluster_QC_stats_sorted_by_rp) %in% c("cluster"))]
  cluster_QC_stats_sorted_by_rp <- cluster_QC_stats_sorted_by_rp[order(cluster_QC_stats_sorted_by_rp$mean_rp, decreasing = TRUE),]
  
  if(find_doublet_info) {
    cluster_QC_stats_sorted_by_doublet <- data.frame("cluster" = cluster_ids, "num_cells" = num_cells, "num_cells_high_viral" = num_cells_high, "num_cells_low_viral" = num_cells_low,
                                                "mean_doublet" = cluster_mean_doublet, "mean_doublet_high_viral" = cluster_mean_doublet_high,
                                                "mean_doublet_low_viral" = cluster_mean_doublet_low, "mean_doublet_diff" = abs(cluster_mean_doublet_high - cluster_mean_doublet_low))
    rownames(cluster_QC_stats_sorted_by_doublet) <- cluster_QC_stats_sorted_by_doublet$cluster
    cluster_QC_stats_sorted_by_doublet <- cluster_QC_stats_sorted_by_doublet[ , !(names(cluster_QC_stats_sorted_by_doublet) %in% c("cluster"))]
    cluster_QC_stats_sorted_by_doublet <- cluster_QC_stats_sorted_by_doublet[order(cluster_QC_stats_sorted_by_doublet$mean_doublet, decreasing = TRUE),]
  }
  if(find_doublet_info) {
    return(list(cluster_distributions, cluster_predictions, cell_cycle_df, cluster_QC_stats_sorted_by_mt, cluster_QC_stats_sorted_by_nFeature, cluster_QC_stats_sorted_by_nCount, cluster_QC_stats_sorted_by_rp, cluster_QC_stats_sorted_by_doublet))
  } else {
    return(list(cluster_distributions, cluster_predictions, cell_cycle_df, cluster_QC_stats_sorted_by_mt, cluster_QC_stats_sorted_by_nFeature, cluster_QC_stats_sorted_by_nCount, cluster_QC_stats_sorted_by_rp))
  }
}

# Create RNA cell type proportion file
create_RNA_cell_type_proportion_file <- function(sc_obj, analysis_dir, group, high_viral_load_samples, d28_samples, male_samples, token = NULL) {
  sample_names <- unique(sc_obj$sample)
  condition_label <- c()
  if(group == "time_point") {
    for(name in sample_names) {
      if(name %in% d28_samples) {
        condition_label <- c(condition_label, "D28")
      } else {
        condition_label <- c(condition_label, "D_minus_1")
      }
    }
  } else if(group == "sex") {
    if(name %in% male_samples) {
      condition_label <- c(condition_label, "MALE")
    } else {
      condition_label <- c(condition_label, "FEMALE")
    }
  } else if(group == "viral_load") {
    for(name in sample_names) {
      if(name %in% high_viral_load_samples) {
        condition_label <- c(condition_label, "HVL")
      } else {
        condition_label <- c(condition_label, "LVL")
      }
    }
  }
  cell_type_proportions_df <- data.frame("Condition" = condition_label, "Sample_name" = paste0("Sample_", sample_names))
  total_cell_counts_df <- data.frame("Sample_name" = paste0("Sample_", sample_names))
  cell_counts <- vector()
  # Find total cell counts for each sample
  for (sample_id in sample_names) {
    idxPass <- which(sc_obj$sample %in% sample_id)
    cellsPass <- names(sc_obj$orig.ident[idxPass])
    sample_subset <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
    cell_counts <- append(cell_counts, ncol(sample_subset))
  }
  total_cell_counts_df <- cbind(total_cell_counts_df, cell_counts)
  for (cell_type in unique(sc_obj$magical_cell_types)) {
    cell_type_proportions <- vector()
    print(cell_type)
    # Grab cells associated with cell type
    idxPass <- which(sc_obj$magical_cell_types %in% cell_type)
    print(length(idxPass))
    cellsPass <- names(sc_obj$orig.ident[idxPass])
    cells_subset <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
    for (sample_id in sample_names) {
      # Subset further based on cells associated with sample ID
      idxPass <- which(cells_subset$sample %in% sample_id)
      cellsPass <- names(cells_subset$orig.ident[idxPass])
      if (length(cellsPass) == 0) {
        cell_type_proportions <- append(cell_type_proportions, 0)
      } else {
        sample_subset <- subset(x = cells_subset, subset = cell_name %in% cellsPass)
        cell_counts <- ncol(sample_subset)
        cell_type_proportions <- append(cell_type_proportions, cell_counts / total_cell_counts_df[total_cell_counts_df$Sample_name == paste0("Sample_", sample_id),]$cell_counts)
      }
    }
    temp_df <- data.frame(cell_type_proportions)
    cell_type <- sub(" ", "_", cell_type)
    names(temp_df)[names(temp_df) == "cell_type_proportions"] <- cell_type
    cell_type_proportions_df <- cbind(cell_type_proportions_df, temp_df)
  }
  if(is.null(token)) {
    write.csv(cell_type_proportions_df, file = paste0(analysis_dir, "RNA_cell_type_proportion_", group, ".csv"), quote = FALSE, row.names = FALSE)  
  } else {
    write.csv(cell_type_proportions_df, file = paste0(analysis_dir, "RNA_cell_type_proportion_", group, "_", token, ".csv"), quote = FALSE, row.names = FALSE)   
  }
}

# Create MAGICAL pseudobulk file
create_magical_cell_type_pseudobulk_files <- function(sc_obj, analysis_dir, token = NULL) {
  print("Computing pseudobulk counts for each cell type")
  for (cell_type in unique(sc_obj$magical_cell_types)) {
    print(cell_type)
    idxPass <- which(sc_obj$magical_cell_types %in% cell_type)
    print(length(idxPass))
    cellsPass <- names(sc_obj$orig.ident[idxPass])
    # Subset again by each sample name
    # Sum by row
    # Keep as vector in list
    cells_subset <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
    cells_pseudobulk <- list()
    for (sample_name in unique(sc_obj$sample)) {
      idxMatch <- which(str_detect(as.character(cells_subset$magical_cell_types),as.character(cell_type)) & str_detect(as.character(cells_subset$sample), sample_name))
      if (length(idxMatch)>1) {
        samples_subset <- subset(x = cells_subset, subset = sample %in% sample_name)
        samples_data <- samples_subset@assays$RNA@counts
        samples_data <- rowSums(as.matrix(samples_data))
        cells_pseudobulk[[sample_name]] <- samples_data
      } else {
        cells_pseudobulk[[sample_name]] <- numeric(nrow(sc_obj@assays$RNA))
      }
    }
    final_cells_pseudobulk_df <- bind_cols(cells_pseudobulk[1])
    for (idx in 2:length(unique(sc_obj$sample))) {
      final_cells_pseudobulk_df <- bind_cols(final_cells_pseudobulk_df, cells_pseudobulk[idx])
    }
    final_cells_pseudobulk_df <- as.data.frame(final_cells_pseudobulk_df)
    rownames(final_cells_pseudobulk_df) <- names(cells_pseudobulk[[1]])
    colnames(final_cells_pseudobulk_df) <- paste0("Sample_", unique(sc_obj$sample))
    cell_type <- sub(" ", "_", cell_type)
    if(is.null(token)) {
      write.csv(final_cells_pseudobulk_df, paste0(analysis_dir, "pseudo_bulk_RNA_count_", cell_type, ".csv"))
    } else { 
      if (!dir.exists(paste0(analysis_dir, token))) {dir.create(paste0(analysis_dir, token), recursive = TRUE)}
      write.csv(final_cells_pseudobulk_df, paste0(analysis_dir, token, "/pseudo_bulk_RNA_count_", cell_type, ".csv"))
    }
  }
}

create_magical_input_files <- function(sc_obj, MAGICAL_file_dir) {
  if (!dir.exists(MAGICAL_file_dir)) {dir.create(MAGICAL_file_dir)}
  MAGICAL_cell_metadata_dir <- paste0(MAGICAL_file_dir, "scRNA_Cell_Metadata/")
  if (!dir.exists(MAGICAL_cell_metadata_dir)) {dir.create(MAGICAL_cell_metadata_dir)}
  MAGICAL_read_counts_dir <- paste0(MAGICAL_file_dir, "scRNA_Read_Counts/")
  if (!dir.exists(MAGICAL_read_counts_dir)) {dir.create(MAGICAL_read_counts_dir)}
  MAGICAL_genes_dir <- paste0(MAGICAL_file_dir, "scRNA_Genes/")
  if (!dir.exists(MAGICAL_genes_dir)) {dir.create(MAGICAL_genes_dir)}
  
  # Write out genes
  write.table(sc_obj@assays$RNA$counts@Dimnames[[1]], file = paste0(MAGICAL_genes_dir, "RNA_genes.tsv"),
              quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")
  
  for(cell_type in unique(sc_obj$magical_cell_types)) {
    print(cell_type)
    cell_type_for_file_name <- sub(" ", "_", cell_type)
    
    #RNA assay cell read counts
    cell_index=which(sc_obj$magical_cell_types==cell_type)
    cell_type_scRNA_counts = as(sc_obj@assays$RNA$counts[, cell_index], "dgTMatrix")
    saveRDS(cell_type_scRNA_counts, file= paste0(MAGICAL_read_counts_dir, cell_type_for_file_name, "_HVL_RNA_read_counts.rds"))
    
    # Cell metadata
    cell_type_scRNA_meta = sc_obj@meta.data[cell_index,]
    write.table(data.frame(rownames(cell_type_scRNA_meta), cell_type_scRNA_meta$magical_cell_types, cell_type_scRNA_meta$sample, cell_type_scRNA_meta$time_point),
                file = paste0(MAGICAL_cell_metadata_dir, cell_type_for_file_name, "_HVL_RNA_cell_metadata.tsv"),
                quote = FALSE, row.names = TRUE, col.names = FALSE,  sep = "\t")
  }
}
