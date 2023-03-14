#------------------------------------------------
# Prototype API for TALOS
#
# Author: Yuan Wang
# Date:  08/08/2022
#------------------------------------------------

SPEEDI_dir <- "~/SPEEDI"
source(paste0(SPEEDI_dir, "/prototype_utils.R"))

# Human CC Genes
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Mouse CC Genes
#cc.gene.updated.mouse <- readRDS(paste0(SPEEDI_dir, "/cc.gene.updated.mouse.rds"))
#m.s.genes <- cc.gene.updated.mouse$m.s.genes
#m.g2m.genes <- cc.gene.updated.mouse$m.g2m.genes

SEED <- 1824409L
set.seed(SEED)

run_SPEEDI <- function(data_path, output_dir, sample_id_list, naming_token, save_progress = FALSE, remove_doublets = FALSE, run_qc = FALSE, run_soup = FALSE) {
  all_sc_exp_matrices <- Read_h5(data_path, sample_id_list)
  assign("sc_obj", FilterRawData(all_sc_exp_matrices, human = TRUE, remove_doublets = FALSE), envir = .GlobalEnv)
  rm(all_sc_exp_matrices)
  # Add cell names as a metadata column - this is handy for selecting subsets of cells
  cell_names <- rownames(sc_obj@meta.data)
  assign("sc_obj", AddMetaData(sc_obj, metadata = cell_names, col.name = "cell_name"), envir = .GlobalEnv)
  # Find intersection between good quality ATAC and our RNA-seq cells
  filtered_ATAC_cells <- read.table("~/high_vs_low_viral_load_D28/ATAC_seq_data_output/ATAC_all_cells.txt", comment.char = "")
  filtered_ATAC_cells <- filtered_ATAC_cells$V1
  sc_obj_ATAC_overlap_all <- subset(sc_obj, cell_name %in% filtered_ATAC_cells)
  assign("sc_obj", InitialProcessing(sc_obj, human = TRUE), envir = .GlobalEnv)
  if(save_progress) {
    save(sc_obj, file = paste0(output_dir, "3_", naming_token, "_sc_obj_full.rds"))
  }
  assign("sc_obj", InferBatches(sc_obj), envir = .GlobalEnv)
  assign("sc_obj", IntegrateByBatch(sc_obj), envir = .GlobalEnv)
  if(save_progress) {
    saveRDS(sc_obj, file = paste0(output_dir, "5_", naming_token, "_sc_obj.rds"))
  }
  assign("sc_obj", VisualizeIntegration(sc_obj), envir = .GlobalEnv)
  if(save_progress) {
    save(sc_obj, file = paste0(output_dir, "6_", naming_token, "_sc_obj.rds"))
  }
  reference <- LoadReference("PBMC", human = TRUE)
  # Remove cell types which we don't believe in
  idx <- which(reference$celltype.l2 %in% c("Doublet", "B intermediate", "CD4 CTL", "gdT", "dnT", "ILC"))
  reference <- reference[,-idx]
  assign("sc_obj", MapCellTypes(sc_obj, reference), envir = .GlobalEnv)
  rm(reference)
  if(save_progress) {
    save(sc_obj, file = paste0(output_dir, "7_", naming_token, "_sc_obj.rds"))
  }
  return(sc_obj)
}

create_seurat_object_with_info <- function(all_sc_exp_matrices) {
  sc_obj <- CreateSeuratObject(counts = all_sc_exp_matrices,
                               assay = "RNA",
                               min.cells = 3,
                               min.features = 3,
                               project = "flu")
  sc_obj$sample <- as.vector(sapply(strsplit(colnames(sc_obj), "#"), "[", 1))
  sc_obj <- PercentageFeatureSet(object = sc_obj,
                                 pattern = "^MT-",
                                 col.name = "percent.mt") 
  sc_obj <- PercentageFeatureSet(object = sc_obj,
                                 pattern = "^RPS",
                                 col.name = "percent.rps") 
  sc_obj <- PercentageFeatureSet(object = sc_obj,
                                 pattern = "^RPL",
                                 col.name = "percent.rpl") 
  sc_obj <- PercentageFeatureSet(object = sc_obj,
                                 pattern = "^HB[A|B]",
                                 col.name = "percent.hb")
  sc_obj <- PercentageFeatureSet(object = sc_obj,
                                 pattern = "^RP[SL]",
                                 col.name = "percent.rp")
  # Adding cell_names as metadata is useful (e.g., for subsetting)
  cell_names <- rownames(sc_obj@meta.data)
  sc_obj <- AddMetaData(sc_obj, metadata = cell_names, col.name = "cell_name")
  return(sc_obj)
}

add_sample_metadata <- function(sc_obj, high_viral_load_samples, low_viral_load_samples,
                                d28_samples, d_minus_1_samples, male_samples, female_samples) {
  # Add viral load metadata
  viral_load_metadata <- vector(length = length(sc_obj$sample))
  idxPass <- which(sc_obj$sample %in% high_viral_load_samples)
  viral_load_metadata[idxPass] <- "HVL"
  idxPass <- which(sc_obj$sample %in% low_viral_load_samples)
  viral_load_metadata[idxPass] <- "LVL"
  sc_obj$viral_load <- viral_load_metadata
  # Add day metadata
  day_metadata <- vector(length = length(sc_obj$sample))
  idxPass <- which(sc_obj$sample %in% d28_samples)
  day_metadata[idxPass] <- "D28"
  idxPass <- which(sc_obj$sample %in% d_minus_1_samples)
  day_metadata[idxPass] <- "D_MINUS_1"
  sc_obj$day <- day_metadata
  # Add sex metadata
  sex_metadata <- vector(length = length(sc_obj$sample))
  idxPass <- which(sc_obj$sample %in% male_samples)
  sex_metadata[idxPass] <- "MALE"
  idxPass <- which(sc_obj$sample %in% female_samples)
  sex_metadata[idxPass] <- "FEMALE"
  sc_obj$sex <- sex_metadata
  # By default, we'll group samples by viral load (high then low)
  sc_obj$sample <- factor(sc_obj$sample, levels = c(high_viral_load_samples, low_viral_load_samples))
  return(sc_obj)
}

process_matrices_through_soup <- function(data_path, sample_id_list) {
  sample_id <- sample_id_list[1]
  # Load data and estimate soup profile
  sc <- load10X(paste0(data_path, sample_id, "/outs/"))
  # Estimate rho
  sc <- autoEstCont(sc, forceAccept = TRUE)
  # Clean the data
  all_sc_exp_matrices <- adjustCounts(sc)
  # Add sample name
  prefix <- paste0(sample_id, "#")
  colnames(all_sc_exp_matrices) <- paste0(prefix, colnames(all_sc_exp_matrices))
  rest_of_sample_list <- sample_id_list[2:length(sample_id_list)]
  for(sample_id in rest_of_sample_list) {
    print(sample_id)
    # Load data and estimate soup profile
    sc = load10X(paste0(data_path, sample_id, "/outs/"))
    # Estimate rho
    sc = autoEstCont(sc, forceAccept = TRUE)
    # Clean the data
    sc_exp_matrix = adjustCounts(sc)
    # Add sample name
    prefix <- paste0(sample_id, "#")
    colnames(sc_exp_matrix) <- paste0(prefix, colnames(sc_exp_matrix))
    # Add to list of all matrices
    all_sc_exp_matrices <- cbind(all_sc_exp_matrices, sc_exp_matrix)
  }
}

generate_qc_plots <- function(all_sc_exp_matrices, plot_dir, date, high_viral_load_samples, low_viral_load_samples,
                                d28_samples, d_minus_1_samples, male_samples, female_samples) {
  sc_obj <- create_seurat_object_with_info(all_sc_exp_matrices)
  sc_obj <- add_sample_metadata(sc_obj, high_viral_load_samples, low_viral_load_samples,
                                d28_samples, d_minus_1_samples, male_samples, female_samples)
  # Set up plotting colors
  viral_load_plotting_colors <- c(rep("#FC4E07", length(high_viral_load_samples)), rep("#2E9FDF", length(low_viral_load_samples)))
  day_plotting_colors <- c(rep("#FC4E07", length(d28_samples)), rep("#2E9FDF", length(d_minus_1_samples)))
  sex_plotting_colors <- c(rep("#FC4E07", length(male_samples)), rep("#2E9FDF", length(female_samples)))
  # Loops for plotting
  group.by.categories <- c("viral_load", "day", "sex")
  plotting.colors <- list(viral_load_plotting_colors, day_plotting_colors, sex_plotting_colors)
  for(index in 1:length(group.by.categories)) {
    current_category <- group.by.categories[index]
    current_plotting_colors <- plotting.colors[[index]]
    # VIRAL LOAD QC PLOTS
    p <- VlnPlot(sc_obj, features = c("nFeature_RNA"), split.by = "sample", group.by = current_category, raster = FALSE) + scale_fill_manual(values = current_plotting_colors) +
      xlab(current_category)
    ggsave(paste0(plot_dir, "nFeature_violin_plots_", date, "_category_", current_category, ".png"), plot = p, device = "png", width = 10, height = 10, units = "in")
    p <- VlnPlot(sc_obj, features = c("nCount_RNA"), split.by = "sample", group.by = current_category, raster = FALSE) + scale_fill_manual(values = current_plotting_colors) +
      xlab(current_category)
    ggsave(paste0(plot_dir, "nCount_violin_plots_", date, "_category_", current_category, ".png"), plot = p, device = "png", width = 10, height = 10, units = "in")
    p <- VlnPlot(sc_obj, features = c("percent.mt"), split.by = "sample", group.by = current_category, raster = FALSE) + scale_fill_manual(values = current_plotting_colors) +
      xlab(current_category)
    ggsave(paste0(plot_dir, "percentMT_violin_plots_", date, "_category_", current_category, ".png"), plot = p, device = "png", width = 10, height = 10, units = "in")
    p <- VlnPlot(sc_obj, features = c("percent.hb"), split.by = "sample", group.by = current_category, raster = FALSE) + scale_fill_manual(values = current_plotting_colors) +
      xlab(current_category)
    ggsave(paste0(plot_dir, "percentHB_violin_plots_", date, "_category_", current_category, ".png"), plot = p, device = "png", width = 10, height = 10, units = "in")
    p <- VlnPlot(sc_obj, features = c("percent.rp"), split.by = "sample", group.by = current_category, raster = FALSE) + scale_fill_manual(values = current_plotting_colors) +
      xlab(current_category)
    ggsave(paste0(plot_dir, "percentRP_violin_plots_", date, "_category_", current_category, ".png"), plot = p, device = "png", width = 10, height = 10, units = "in")
  }

  # Create nCount vs nFeature scatter plot for each sample
  individual_samples <- SplitObject(sc_obj, split.by = "sample")
  # Visualize nCount vs nFeature via scatter plot for each sample
  nCount_vs_nFeature_plots <- list()
  for (i in 1:length(individual_samples)) {
    p <- FeatureScatter(individual_samples[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")
    nCount_vs_nFeature_plots[[i]] <- p
  }
  n <- length(nCount_vs_nFeature_plots)
  nCol <- floor(sqrt(n))
  nCount_vs_nFeature_plots <- do.call("grid.arrange", c(nCount_vs_nFeature_plots, ncol=nCol))
  ggsave(paste0(plot_dir, "nCount_vs_nFeature_plots_", date,".png"), plot = nCount_vs_nFeature_plots, device = "png", width = 20, height = 20, units = "in")
  rm(individual_samples)
  rm(sc_obj)
}

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

combine_cell_types_magical <- function(sc_obj, resolution = 1.5) {
  Cell_type_combined <- sc_obj$predicted_celltype_majority_vote
  levels(Cell_type_combined) <- c(levels(Cell_type_combined), "T Naive", "B")
  idx <- grep("CD4 Naive", Cell_type_combined)
  Cell_type_combined[idx] <- "T Naive"
  idx <- grep("CD8 Naive", Cell_type_combined)
  Cell_type_combined[idx] <- "T Naive"
  idx <- grep("Treg", Cell_type_combined)
  Cell_type_combined[idx] <- "T Naive"
  idx <- grep("NK_CD56bright", Cell_type_combined)
  Cell_type_combined[idx] <- "NK"
  idx <- grep("B", Cell_type_combined)
  Cell_type_combined[idx] <- "B"
  sc_obj$magical_cell_types <- Cell_type_combined
  return(sc_obj)
}

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
  all_cell_types <- union(unique(sc_obj$predicted_celltype_majority_vote), unique(sc_obj$magical_cell_types))
  for (cell_type in all_cell_types) {
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
    print(table(cells_subset$viral_load))
    DefaultAssay(cells_subset) <- "SCT"
    Idents(cells_subset) <- group
    if(group == "viral_load") {
      first_group <- "HIGH"
      second_group <- "LOW"
    } else if(group == "day") {
      first_group <- "D28"
      second_group <- "D_minus_1"
    } else if(group == "sex") {
      first_group <- "MALE"
      second_group <- "FEMALE"
    }
    diff_markers <- FindMarkers(cells_subset, ident.1 = first_group, ident.2 = second_group, assay = "SCT", recorrect_umi = FALSE, logfc.threshold = 0, min.pct = 0)
    cell_type <- sub(" ", "_", cell_type)
    write.csv(diff_markers, paste0(analysis_dir, "HIGH-vs-LOW-degs-", cell_type, ".csv"), quote = FALSE)
  }
}

create_magical_cell_type_proportion_file <- function(sc_obj, group) {
  cell_type_proportions_df <- data.frame("Condition" = viral_load_label, "Sample_name" = all_viral_load)
  total_cell_counts_df <- data.frame("Sample_name" = all_viral_load)
  cell_counts <- vector()
  # Find total cell counts for each sample
  for (sample_id in all_viral_load) {
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
    for (sample_id in all_viral_load) {
      # Subset further based on cells associated with sample ID
      idxPass <- which(cells_subset$sample %in% sample_id)
      cellsPass <- names(cells_subset$orig.ident[idxPass])
      if (length(cellsPass) == 0) {
        cell_type_proportions <- append(cell_type_proportions, 0)
      } else {
        sample_subset <- subset(x = cells_subset, subset = cell_name %in% cellsPass)
        cell_counts <- ncol(sample_subset)
        cell_type_proportions <- append(cell_type_proportions, cell_counts / total_cell_counts_df[total_cell_counts_df$Sample_name == sample_id,]$cell_counts)
      }
    }
    temp_df <- data.frame(cell_type_proportions)
    names(temp_df)[names(temp_df) == "cell_type_proportions"] <- cell_type
    cell_type_proportions_df <- cbind(cell_type_proportions_df, temp_df)
  }
  write.csv(cell_type_proportions_df, file = paste0(output_dir, "RNA_cell_type_proportion.csv"), quote = FALSE, row.names = FALSE)
}

print_UMAP <- function(sc_obj, sample_count, group_by_category, plot_dir, file_name) {
  cell_count <- length(sc_obj$cell_name)
  current_title <- paste0("snRNA-seq Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
  DimPlot(sc_obj, reduction = "umap", group.by = group_by_category, label = TRUE,
        label.size = 3, repel = TRUE, raster = FALSE) +
  labs(title = current_title) +
  theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(plot_dir, file_name), device = "png", dpi = 300)
}

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

capture_cluster_info <- function(sc_obj) {
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
    #cluster_mean_doublet <- c(cluster_mean_doublet, mean(filtered_cluster$scDblFinder.score))
    cluster_prediction <- sort(table(filtered_cluster$predicted_celltype_majority_vote), decreasing = TRUE)[1]
    cluster_predictions <- append(cluster_predictions, cluster_prediction)
    cluster_distributions[[idx]] <- table(filtered_cluster$predicted.id)
    cluster_mean_S_score <- append(cluster_mean_S_score, mean(filtered_cluster$S.Score))
    cluster_mean_G2M_score <- append(cluster_mean_G2M_score, mean(filtered_cluster$G2M.Score))
    cluster_mean_CC_difference <- append(cluster_mean_CC_difference, mean(filtered_cluster$CC.Difference))
    # High
    idxPass <- which(filtered_cluster$viral_load %in% "HVL")
    cellsPass <- names(filtered_cluster$orig.ident[idxPass])
    num_cells_high <- c(num_cells_high, length(cellsPass))
    if(length(cellsPass) > 0) {
      filtered_cluster_high <- subset(x = filtered_cluster, subset = cell_name %in% cellsPass)
      cluster_mean_mito_high <- c(cluster_mean_mito_high, mean(filtered_cluster_high$percent.mt))
      cluster_mean_nFeature_high <- c(cluster_mean_nFeature_high, mean(filtered_cluster_high$nFeature_RNA))
      cluster_mean_nCount_high <- c(cluster_mean_nCount_high, mean(filtered_cluster_high$nCount_RNA))
      cluster_mean_rp_high <- c(cluster_mean_rp_high, mean(filtered_cluster_high$percent.rp))
      #cluster_mean_doublet_high <- c(cluster_mean_doublet_high, mean(filtered_cluster_high$scDblFinder.score))
    } else {
      cluster_mean_mito_high <- c(cluster_mean_mito_high, NA)
      cluster_mean_nFeature_high <- c(cluster_mean_nFeature_high, NA)
      cluster_mean_nCount_high <- c(cluster_mean_nCount_high, NA)
      cluster_mean_rp_high <- c(cluster_mean_rp_high, NA)
      #cluster_mean_doublet_high <- c(cluster_mean_doublet_high, NA)
    }
    # Low
    idxPass <- which(filtered_cluster$viral_load %in% "LVL")
    cellsPass <- names(filtered_cluster$orig.ident[idxPass])
    num_cells_low <- c(num_cells_low, length(cellsPass))
    if(length(cellsPass) > 0) {
      filtered_cluster_low <- subset(x = filtered_cluster, subset = cell_name %in% cellsPass)
      cluster_mean_mito_low <- c(cluster_mean_mito_low, mean(filtered_cluster_low$percent.mt))
      cluster_mean_nFeature_low <- c(cluster_mean_nFeature_low, mean(filtered_cluster_low$nFeature_RNA))
      cluster_mean_nCount_low <- c(cluster_mean_nCount_low, mean(filtered_cluster_low$nCount_RNA))
      cluster_mean_rp_low <- c(cluster_mean_rp_low, mean(filtered_cluster_low$percent.rp))
      #cluster_mean_doublet_low <- c(cluster_mean_doublet_low, mean(filtered_cluster_low$scDblFinder.score))
    } else {
      cluster_mean_mito_low <- c(cluster_mean_mito_low, NA)
      cluster_mean_nFeature_low <- c(cluster_mean_nFeature_low, NA)
      cluster_mean_nCount_low <- c(cluster_mean_nCount_low, NA)
      cluster_mean_rp_low <- c(cluster_mean_rp_low, NA)
      #cluster_mean_doublet_low <- c(cluster_mean_doublet_low, NA)
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
  
  
  #cluster_QC_stats_sorted_by_doublet <- data.frame("cluster" = cluster_ids, "num_cells" = num_cells, "num_cells_high_viral" = num_cells_high, "num_cells_low_viral" = num_cells_low,
  #                                            "mean_doublet" = cluster_mean_doublet, "mean_doublet_high_viral" = cluster_mean_doublet_high,
  #                                            "mean_doublet_low_viral" = cluster_mean_doublet_low, "mean_doublet_diff" = abs(cluster_mean_doublet_high - cluster_mean_doublet_low))
  #rownames(cluster_QC_stats_sorted_by_doublet) <- cluster_QC_stats_sorted_by_doublet$cluster
  #cluster_QC_stats_sorted_by_doublet <- cluster_QC_stats_sorted_by_doublet[ , !(names(cluster_QC_stats_sorted_by_doublet) %in% c("cluster"))]
  #cluster_QC_stats_sorted_by_doublet <- cluster_QC_stats_sorted_by_doublet[order(cluster_QC_stats_sorted_by_doublet$mean_doublet, decreasing = TRUE),]
  return(list(cluster_distributions, cluster_predictions, cell_cycle_df, cluster_QC_stats_sorted_by_mt, cluster_QC_stats_sorted_by_nFeature, cluster_QC_stats_sorted_by_nCount, cluster_QC_stats_sorted_by_rp))
}


Read_h5 <- function(data_path, sample_id_list) {
  message("Step 1: Reading all samples...")

  # Make reading data parallel
  if (Sys.getenv("SLURM_NTASKS_PER_NODE") == "") {
    n.cores <- as.numeric(detectCores())
  } else {
    n.cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
  }

  if (n.cores > length(sample_id_list)) {
    n.cores <- length(sample_id_list)
  }

  message(paste0("Number of cores: ", n.cores))

  registerDoMC(n.cores)
  message("Begin parallelizing...")

  all_sc_exp_matrices <- foreach(
    i = 1:length(sample_id_list),
    .combine = 'cbind',
    .packages = c("Seurat", "base")
  ) %dopar% {
    library(hdf5r)
#     print(paste0(sample_id_list[[i]], "/filtered_feature_bc_matrix"))
#     sc_matrix <- Read10X(paste0(sample_id_list[[i]], "/filtered_feature_bc_matrix"))

    print(paste0(data_path, sample_id_list[[i]], "/outs/filtered_feature_bc_matrix.h5"))
    sc_matrix <- Read10X_h5(paste0(data_path, sample_id_list[[i]], "/outs/filtered_feature_bc_matrix.h5"))

    if (class(x = sc_matrix) == "list") {
      sc_exp_matrix <- sc_matrix$`Gene Expression`
    } else {
      sc_exp_matrix <- sc_matrix
    }
    if (grepl("_|\\.", i)) {
      prefix <- paste0(strsplit(sample_id_list[[i]], "_")[[1]][1], "#")
    } else {
      prefix <- paste0(sample_id_list[[i]], "#")
    }
    colnames(sc_exp_matrix) <- paste0(prefix, colnames(sc_exp_matrix))
    return(sc_exp_matrix)
  }

  message(paste0("Raw data has ", dim(all_sc_exp_matrices)[2], " barcodes and ", dim(all_sc_exp_matrices)[1], " transcripts."))
  return(all_sc_exp_matrices)
}

FilterRawData <- function(sc_obj, human = TRUE, record_doublets = FALSE, adaptive_QC_thresholds = TRUE, 
                          min_nFeature = 900, max_nFeature = 4000, max_nCount = 10000, max_percent_mt = 15, 
                          max_percent_hb = 0.4, max_percent_rp = 50) {
  message("Step 2: Filtering out bad samples...")

  sc_obj <- CreateSeuratObject(counts = all_sc_exp_matrices,
                               assay = "RNA",
                               min.cells = 3,
                               min.features = 3,
                               project = "flu")

  sc_obj$sample <- as.vector(sapply(strsplit(colnames(sc_obj), "#"), "[", 1))
  # Adding cell_names as metadata is useful (e.g., for subsetting)
  cell_names <- rownames(sc_obj@meta.data)
  sc_obj <- AddMetaData(sc_obj, metadata = cell_names, col.name = "cell_name")

  if(record_doublets) {
    message("Recording doublets...")
    sc_obj <- Seurat::as.Seurat(scDblFinder::scDblFinder(Seurat::as.SingleCellExperiment(sc_obj), samples = "sample", BPPARAM=MulticoreParam(7, RNGseed=SEED)))
    # See distribution of doublets in each sample
    doublet_sc_obj <- subset(x = sc_obj, subset = scDblFinder.class %in% "doublet")
    message("Number of doublets removed in each sample:")
    print(table(doublet_sc_obj$sample))
    rm(doublet_sc_obj)
  }

  if (human) {
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^MT-",
                                   col.name = "percent.mt")
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^RPS",
                                   col.name = "percent.rps")
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^RPL",
                                   col.name = "percent.rpl")
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^HB[A|B]",
                                   col.name = "percent.hb")
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^RP[SL]",
                                   col.name = "percent.rp")


  } else {
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^mt-",
                                   col.name = "percent.mt")
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^Rps",
                                   col.name = "percent.rps")
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^Rpl",
                                   col.name = "percent.rpl")
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^Hb[a|b]",
                                   col.name = "percent.hb")
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^Rp[sl]",
                                   col.name = "percent.rp")
  }

  objects <- SplitObject(sc_obj, split.by = "sample")


  if (Sys.getenv("SLURM_NTASKS_PER_NODE") == "") {
    n.cores <- as.numeric(detectCores())
  } else {
    n.cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
  }

  if (n.cores > length(objects)) {
    n.cores <- length(objects)
  }

  message(paste0("Number of cores: ", n.cores))

  registerDoMC(n.cores)
  message("Begin parallizing...")

  sc_obj <- foreach(
    i = 1:length(objects),
    .combine = 'merge',
    .packages = c("Seurat", "base")
  ) %dopar% {

#     lower_nF_dec <- kneedle(hist(objects[[i]]$nFeature_RNA, breaks=100, plot=F)$breaks[-1],
#                             hist(objects[[i]]$nFeature_RNA, breaks=100, plot=F)$counts,
#                             decreasing = F)[1]
#     lower_nF_inc <- kneedle(hist(objects[[i]]$nFeature_RNA, breaks=100, plot=F)$breaks[-1],
#                             hist(objects[[i]]$nFeature_RNA, breaks=100, plot=F)$counts,
#                             decreasing = T)[1]
#     lower_nF <- min(lower_nF_dec, lower_nF_inc)
    current_sample_name <- unique(objects[[i]]$sample)
    if(adaptive_QC_thresholds) {
      lower_nF <- kneedle(hist(objects[[i]]$nFeature_RNA, breaks=100, plot=F)$breaks[-1],
                        hist(objects[[i]]$nFeature_RNA, breaks=100, plot=F)$counts)[1]
      if (lower_nF > 1000) { lower_nF <- 1000 }

      if (max(objects[[i]]$percent.mt) > 0) {
         if (max(objects[[i]]$percent.mt) < 5) {
              max_mt <- quantile(objects[[i]]$percent.mt, .99)
          } else {
           max_mt <- kneedle(hist(objects[[i]]$percent.mt, breaks=max(10, 0.5 * max(objects[[i]]$percent.mt)), plot=F)$breaks[-1],
                              hist(objects[[i]]$percent.mt, breaks=max(10, 0.5 * max(objects[[i]]$percent.mt)), plot=F)$counts)[1]
            max_mt <- max(max_mt, quantile(objects[[i]]$percent.mt, .75))
          }
     } else { max_mt <- 0}

     max_hb <- quantile(objects[[i]]$percent.hb, .99)
     if (max_hb > 10) { max_hb <- 10 }

     object <- subset(x = objects[[i]],
                       subset = nFeature_RNA >= lower_nF &
                         nFeature_RNA < quantile(objects[[i]]$nFeature_RNA, .99) &
                        percent.mt <= max_mt &
                         percent.rps <= quantile(objects[[i]]$percent.rps, .99) &
                        percent.rpl <= quantile(objects[[i]]$percent.rpl, .99) &
                         percent.hb <= max_hb)
    } else {
      object <- subset(objects[[i]], nFeature_RNA > min_nFeature & nFeature_RNA < max_nFeature & nCount_RNA < max_nCount & percent.mt < max_percent_mt & percent.hb < max_percent_hb & percent.rp < max_percent_rp)
    }
    if(adaptive_QC_thresholds) {
      print(paste0("THRESHOLDS USED FOR SAMPLE", current_sample_name))
      print(paste0("lower nFeature: ", lower_nF))
      print(paste0("upper nFeature: ", quantile(objects[[i]]$nFeature_RNA, .99)))
      print(paste0("max mt: ", max_mt))
      print(paste0("max rps: ", quantile(objects[[i]]$percent.rps, .99)))
      print(paste0("max rpl: ", quantile(objects[[i]]$percent.rpl, .99)))
      print(paste0("max hb: ", max_hb))
    }

    return(object)
  }

  message(paste0("Filtered data has ", dim(sc_obj)[2], " barcodes and ", dim(sc_obj)[1], " transcripts."))
  return(sc_obj)
}

InitialProcessing <- function(sc_obj, human) {
  message("Step 3: Processing raw data...")
  if (human) {
      sc_obj <- CellCycleScoring(object = sc_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  } else {
       sc_obj <- CellCycleScoring(object = sc_obj, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
  }
  sc_obj$CC.Difference <- sc_obj$S.Score - sc_obj$G2M.Score
  system.time(sc_obj <- SCTransform(object = sc_obj,
                                    vst.flavor = "v2",
                                    vars.to.regress = c("percent.mt",
                                                        "percent.rps",
                                                        "percent.rpl",
                                                        "percent.hb",
                                                        "CC.Difference"),
                                    do.scale = TRUE,
                                    do.center = TRUE,
                                    return.only.var.genes = TRUE,
                                    seed.use = SEED,
                                    verbose = TRUE))
#   sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 3000)
#   sc_obj <- ScaleData(sc_obj)
  sc_obj <- RunPCA(sc_obj, npcs = 30, approx = T, verbose = T, seed.use = SEED)
  sc_obj <- RunUMAP(sc_obj, reduction = "pca", dims = 1:30, seed.use = SEED)
  return(sc_obj)
}

InferBatches <- function(sc_obj) {
  message("Step 4: Infer heterogeneous groups for integration...")

  sc_obj <- FindNeighbors(object = sc_obj, dims = 1:30)
  sc_obj <- FindClusters(object = sc_obj, resolution = 0.1, algorithm = 2, random.seed = SEED)
#   if (length(levels(sc_obj$seurat_clusters)) > 20) {
#       sc_obj <- FindClusters(object = sc_obj, resolution = 0.02, algorithm = 2, random.seed = SEED)
#   }


  # Use LISI metric to guess batch labels
  X <- sc_obj@reductions$umap@cell.embeddings
  meta_data <- data.frame(sc_obj$sample)
  colnames(meta_data) <- "batch"
  meta_data$cluster <- sc_obj$seurat_clusters
  lisi.res <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(lisi.res) <- c("batch", "score", "cluster", "freq")
  clusters.interest <- names(table(sc_obj$seurat_clusters))[prop.table(table(sc_obj$seurat_clusters)) > 0.01]
  for (cluster in clusters.interest) { #levels(sc_obj$seurat_clusters)) {
    cells <- names(sc_obj$seurat_clusters[sc_obj$seurat_clusters == cluster])
    X.sub <- X[which(rownames(X) %in% cells),]
    meta_data.sub <- meta_data[which(rownames(meta_data) %in% cells),]
    res <- compute_lisi(X.sub, meta_data.sub, label_colnames = "batch")
    rownames(res) <- cells
    colnames(res) <- "score"
    res$batch <- meta_data.sub$batch
    agg.res <- aggregate(.~batch,data=res,mean)
    agg.res$cluster <- cluster
    agg.res$freq <- data.frame(table(res$batch))$Freq[which(data.frame(table(res$batch))$Var1 %in% agg.res$batch)]
    lisi.res <- rbind(lisi.res, agg.res)
  }

  p.values <- list()
  used.sample.dump <- c()
  batch.assign <- list()
  for ( i in clusters.interest) {
    #print(i)
    #lisi.res.sub <- lisi.res[lisi.res$score <= quantile(lisi.res$score, 1),]
    lisi.res.sub <- lisi.res[lisi.res$cluster == i,]
    if (max(lisi.res.sub$score) <= 1.1) {
      samples.of.batch <- lisi.res.sub$batch[1]
      if (!(samples.of.batch %in% used.sample.dump)) {
        batch.assign <- lappend(batch.assign, samples.of.batch)
      }
      used.sample.dump <- union(used.sample.dump, samples.of.batch)
    } else {
      lisi.res.sub$scaled.score <- scale_zero_one(lisi.res.sub$score * (lisi.res.sub$freq / sum(lisi.res.sub$freq)))
      lisi.res.sub <- lisi.res.sub[order(lisi.res.sub$scaled.score, decreasing = TRUE),]
      if (dim(lisi.res.sub)[1] > 30) {
        lisi.res.sub <- lisi.res.sub[1:30,]
      }
      lisi.res.sub$diff.scaled.score <- abs(c(diff(lisi.res.sub$scaled.score), 0))

      if (dim(lisi.res.sub)[1] >= 3) {
        p.values[[i]] <- dixon.test(lisi.res.sub$diff.scaled.score)$p.value[[1]]
      } else {
        p.values[[i]] <- 1
      }

      if (p.values[[i]] < 0.05) {
        max.index <- which.max(lisi.res.sub$diff.scaled.score)
        samples.of.batch <- lisi.res.sub$batch[1:max.index]

        if (any(samples.of.batch %in% used.sample.dump)) {
          if (!all(samples.of.batch %in% used.sample.dump)) {
            used.index <- which(samples.of.batch %in% used.sample.dump)
            samples.of.batch <- samples.of.batch[-used.index]
            if (length(samples.of.batch) > 0) {
                batch.assign <- lappend(batch.assign, samples.of.batch)
            }
          } else if (!list(samples.of.batch) %in% batch.assign) {
              if (length(samples.of.batch) == 1) {
                  batch.assign <- lappend(batch.assign, samples.of.batch)
              } else {
                  used.index <- which(samples.of.batch %in% unlist(batch.assign))
                  samples.of.batch <- samples.of.batch[-used.index]
                  if (length(samples.of.batch) > 0) {
                      batch.assign <- lappend(batch.assign, samples.of.batch)
                  }
              }
          }
        } else {
          batch.assign <- lappend(batch.assign, samples.of.batch)
        }
        used.sample.dump <- union(used.sample.dump, samples.of.batch)
      }
    }
  }

  batch <- as.factor(sc_obj$sample)

  if (length(batch.assign) > 0) {
      levels.batch <- levels(batch)
      for (i in 1:length(batch.assign)) {
          levels.batch[which(levels(batch) %in% batch.assign[[i]])] <- i
      }
      levels.batch[!levels.batch %in% c(1:length(batch.assign))] <- length(batch.assign)+1
      levels(batch) <- levels.batch
      sc_obj$batch <- as.character(batch)
  }

  else {
      message("No batch effect detected!")
      sc_obj$batch <- "No Batch"
  }

  print(unique(batch))

#   saveRDS(sc_obj, paste0(SPEEDI_dir, "/unintegrated.object.RData"))
  return(sc_obj)
}

IntegrateByBatch <- function(sc_obj) {
  message("Step 5: Integrate samples based on inferred groups...")
  sc_obj_list <- SplitObject(sc_obj, split.by = "batch")


  if (Sys.getenv("SLURM_NTASKS_PER_NODE") == "") {
    n.cores <- as.numeric(detectCores())
  } else {
    n.cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
  }

  if (n.cores > length(sc_obj_list)) {
    n.cores <- length(sc_obj_list)
  }

  message(paste0("Number of cores: ", n.cores))

  registerDoMC(n.cores)
  message("Begin parallizing...")


  r <- foreach(
    i = 1:length(sc_obj_list),
    .combine = 'c',
    .packages = c("Seurat", "base")
  ) %dopar% {
    SEED <- 1824409L
      tmp <- SCTransform(object = sc_obj_list[[i]],
                       vst.flavor = "v2",
                       vars.to.regress = c("percent.mt",
                                           "percent.rps",
                                           "percent.rpl",
                                           "percent.hb",
                                           "CC.Difference"),
                       do.scale = TRUE,
                       do.center = TRUE,
                       return.only.var.genes = TRUE,
                       seed.use = SEED,
                       verbose = TRUE)
      tmp <- RunPCA(tmp, npcs = 30, approx = T, verbose = T, seed.use = SEED)
      return(tmp)
  }
  message(paste0(length(r), " samples transformed."))
  message("... Done parallizing")


  message("Select integration features...")
  features <- SelectIntegrationFeatures(object.list = r, nfeatures = 3000)
  r <- PrepSCTIntegration(object.list = r, anchor.features = features)

  message("Find integration anchors...")


#   anchors <- FindIntegrationAnchors(object.list = r,
#                                     normalization.method = "SCT",
#                                     anchor.features = features)

#   if (length(sc_obj_list) > 10) {
#     message("...use reference-based integration...")
#     anchors <- FindIntegrationAnchors(object.list = r,
#                                       reference = 1,
#                                       normalization.method = "SCT",
#                                       anchor.features = features,
#                                       reduction = "rpca",
#                                       k.anchor = 10)
#   } else {
    anchors <- FindIntegrationAnchors(object.list = r,
                                      normalization.method = "SCT",
                                      anchor.features = features,
                                      reduction = "rpca",
                                      k.anchor = 10)
#    }

  message("Begin integration...")
#   integrated_obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  integrated_obj <- IntegrateData(anchorset = anchors,
                                  normalization.method = "SCT",
                                  k.weight = 100)
  DefaultAssay(integrated_obj) <- "integrated"

  rm(sc_obj_list)
  rm(features)
  rm(anchors)

  return(integrated_obj)
  # return(r)
}




VisualizeIntegration <- function(sc_obj, prep_sct_find_markers = TRUE) {
  set.seed(1824409L)
  sc_obj <- ScaleData(sc_obj, verbose = T)
#   sc_obj <- PCA(sc_obj)
  sc_obj <- RunPCA(sc_obj, npcs = 30, approx = T, verbose = T, seed.use = SEED)
  sc_obj <- RunUMAP(sc_obj, reduction = "pca", dims = 1:30, seed.use = SEED, return.model = T)
#   sc_obj <- RunUMAP(sc_obj, reduction = "prcomp", dims = 1:30, seed.use = SEED, return.model = T)
  DefaultAssay(sc_obj) <- "SCT"
  if(prep_sct_find_markers) {
    sc_obj <- PrepSCTFindMarkers(sc_obj)
  }
  return(sc_obj)
}

LoadReference <- function(tissue, human) {
  if (human) {
    if (tissue == "Adipose") {
        InstallData("adiposeref")
        return(data("adiposeref")) }

    if (tissue == "Bone Marrow") {
        InstallData("bonemarrowref")
        return(data("bonemarrowref")) }

    if (tissue == "Fetus") {
        InstallData("fetusref")
        return(data("fetusref")) }

    if (tissue == "Heart") {
        InstallData("heartref")
        return(data("heartref")) }

    if (tissue == "Cortex") {
        InstallData("humancortexref")
        return(data("humancortexref")) }

    if (tissue == "Kidney") {
        InstallData("kidneyref")
        return(data("kidneyref")) }

    if (tissue == "Lung") {
        InstallData("lungref")
        return(data("lungref")) }

    if (tissue == "Pancreas") {
        InstallData("pancreasref")
        return(data("pancreasref")) }

    if (tissue == "PBMC") {
        reference <- LoadH5Seurat(paste0("~/SPEEDI/reference/pbmc_multimodal.h5seurat"))
        return(reference) }

    if (tissue == "Tonsil") {
        InstallData("tonsilref")
        return(data("tonsilref")) }
  }
  if (!human) {
    if (tissue == "Cortex") {
        InstallData("mousecortexref")
        return(data("mousecortexref")) }
  }
}

FindMappingAnchors <- function(sc_obj, reference, data_type = "scRNA") {
  DefaultAssay(sc_obj) <- "integrated"
  if(data_type == "scRNA") {
    recompute.residuals.value <- "T"
  } else {
    recompute.residuals.value <- "F"
  }
  anchors <- FindTransferAnchors(reference = reference,
                                 query = sc_obj,
                                 normalization.method = "SCT",
                                 recompute.residuals = recompute.residuals.value,
                                 reference.reduction = "spca")
  return(anchors)
}

MajorityVote <- function(sc_obj, current_resolution = 1.5) {
  associated_res_attribute <- paste0("integrated_snn_res.", current_resolution)
  message("Begin majority voting...")
  DefaultAssay(sc_obj) <- "integrated"
  sc_obj <- FindNeighbors(sc_obj, reduction = "pca", dims = 1:30)
  # TODO: Add code to find the best resolution (e.g., by using Clustree?)
  sc_obj <- FindClusters(sc_obj, resolution = current_resolution)
  sc_obj$predicted.id <- as.character(sc_obj$predicted.id)

  #idx <- grep("CD4 T", sc_obj$predicted.id)
  #sc_obj$predicted.id[idx] <- "CD4 Memory"
  #idx <- grep("CD8 T", sc_obj$predicted.id)
  #sc_obj$predicted.id[idx] <- "CD8 Memory"
  #idx <- grep("cDC", sc_obj$predicted.id)
  #sc_obj$predicted.id[idx] <- "cDC"
  #idx <- grep("Proliferating", sc_obj$predicted.id)
  #sc_obj$predicted.id[idx] <- "Proliferating"
  #idx <- grep("B", sc_obj$predicted.id)
  #sc_obj$predicted.id[idx] <- "B"
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "CD4 Naive", "T Naive")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "CD8 Naive", "T Naive")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "NK_CD56bright", "NK")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "ASDC", "CD14 Mono")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "cDC", "CD14 Mono")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "Eryth", "CD14 Mono")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "HSPC", "CD14 Mono")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "pDC", "CD14 Mono")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "Plasmablast", "CD14 Mono")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "Platelet", "CD14 Mono")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "Treg", "T Naive")

  integrated_snn_res_df <- sc_obj[[associated_res_attribute]]
  integrated_snn_res_cell_names <- rownames(integrated_snn_res_df)
  integrated_snn_res_values <- integrated_snn_res_df[,1]

  cluster.dump <- as.numeric(levels(integrated_snn_res_values))
  sc_obj$predicted_celltype_majority_vote <- sc_obj$seurat_clusters
  levels(sc_obj$predicted_celltype_majority_vote) <- as.character(levels(sc_obj$predicted_celltype_majority_vote))
  for (i in unique(sc_obj$predicted.id)) {
    print(i)
    cells <- names(sc_obj$predicted.id[sc_obj$predicted.id == i])
    freq.table <- as.data.frame(table(integrated_snn_res_df[cells,]))
    freq.table <- freq.table[order(freq.table$Freq, decreasing = TRUE),]
    freq.table$diff <- abs(c(diff(freq.table$Freq), 0))
    if(nrow(freq.table) > 30) {
      freq.table <- freq.table[1:30,]
    }
    p.values <- dixon.test(freq.table$diff)$p.value[[1]]
    max.index <- which.max(freq.table$diff)
    clusters <- as.numeric(as.character(freq.table$Var1[1:max.index]))
    levels(sc_obj$predicted_celltype_majority_vote)[levels(sc_obj$predicted_celltype_majority_vote) %in% as.character(clusters)] <- i
    cluster.dump <- cluster.dump[!cluster.dump %in% clusters]
  }

  if (length(cluster.dump) > 0) {
      for (i in cluster.dump) {
          cells <- rownames(subset(integrated_snn_res_df, integrated_snn_res_df[,1] == i,))
          freq.table <- as.data.frame(table(sc_obj$predicted.id[cells]))
          levels(sc_obj$predicted_celltype_majority_vote)[levels(sc_obj$predicted_celltype_majority_vote) %in% as.character(i)] <- as.vector(freq.table$Var1)[which.max(freq.table$Freq)]
      }
  }


  message("...End majority voting")
  return(sc_obj)
}

MapCellTypes <- function(sc_obj, reference, data_type = "scRNA") {
  message("Step 6: Reference-based cell type mapping...")
  anchors <- FindMappingAnchors(sc_obj, reference, data_type)
  sc_obj <- MapQuery(anchorset = anchors,
                     query = sc_obj,
                     reference = reference,
                     refdata = "celltype.l2",
                     reference.reduction = "spca",
                     reduction.model = "wnn.umap",
                     verbose = TRUE)
  sc_obj <- MajorityVote(sc_obj)
  return(sc_obj)
}




