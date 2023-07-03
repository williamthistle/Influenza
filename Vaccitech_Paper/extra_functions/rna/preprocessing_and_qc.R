normalize_dir_path <- function(path) {
  path <- normalizePath(path, "/")
  # Add "/" to end of path if not already present
  last_char_of_path <- substr(path, nchar(path), nchar(path))
  if(last_char_of_path != "/") {
    path <- paste0(path, "/")
  }
  return(path)
}

print_SPEEDI <- function(current_message, log_flag = FALSE, silence_time = FALSE) {
  if(!silence_time) {
    current_message <- paste0(Sys.time(), ": ", current_message)
  }
  message(current_message)
  if(log_flag) {
    logr::log_print(current_message, console = FALSE, hide_notes = TRUE, blank_after = FALSE)
  }
  return(TRUE)
}

# Create Seurat object with important info for QC
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

# Add metadata for each sample to Seurat obj
add_sample_metadata <- function(sc_obj, high_viral_load_samples, low_viral_load_samples,
                                d28_samples, d_minus_1_samples, male_samples, female_samples, sample_metadata) {
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
  # Add subject metadata
  subject_metadata <- vector(length = length(sc_obj$sample))
  for(current_row in 1:nrow(sample_metadata)) {
    current_sample <- sample_metadata[current_row,]$aliquot
    current_subject <- sample_metadata[current_row,]$subject_id
    idxPass <- which(sc_obj$sample %in% current_sample)
    subject_metadata[idxPass] <- current_subject
  }
  sc_obj$subject <- subject_metadata
  # By default, we'll group samples by viral load (high then low)
  sc_obj$sample <- factor(sc_obj$sample, levels = c(high_viral_load_samples, low_viral_load_samples))
  return(sc_obj)
}

# Process sc/snRNA-seq count matrices through soupX to remove ambient RNA
# TODO: Make this parallel
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

# Generate QC plots for sc/snRNA-seq data
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

# Process all samples individually to get cleaner plots of each one (and make sure integration didn't mess anything up)
process_samples_individually <- function(data_path, sample_id_list, reference, analysis_dir) {
  n.cores <- 8
  registerDoMC(n.cores)
  sc_obj <- foreach(
    i = 1:length(sample_id_list),
    .combine = 'cbind',
    .packages = c("Seurat", "base")
  ) %dopar% {
    current_sample <- sample_id_list[[i]]
    current_matrix <- Read_h5(data_path, current_sample)
    current_obj <- FilterRawData(current_matrix, human = TRUE, record_doublets = FALSE, adaptive_QC_thresholds = FALSE)
    sample_name <- unique(current_obj$sample)
    print(sample_name)
    current_obj <- InitialProcessing(current_obj, human = TRUE)
    current_obj <- MapCellTypes(current_obj, reference, data_type = "snRNA")
    current_obj <- combine_cell_types_initial(current_obj)
    print_UMAP(current_obj, 1, "predicted_celltype_majority_vote", plot_dir, paste0("pre_", sample_name, "_majority_vote_", date, ".png"))
    print_UMAP(current_obj, 1, "seurat_clusters", plot_dir, paste0("pre_", sample_name, "_clusters_", date, ".png"))
    print_UMAP(current_obj, 1, "predicted.id", plot_dir, paste0("pre_", sample_name, "_raw_predictions_", date, ".png"))
    save(current_obj, file = paste0(analysis_dir, paste0(sample_name, ".rds")))
  }
}

#process_samples_individually(data_path, sample_id_list, reference, analysis_dir)
