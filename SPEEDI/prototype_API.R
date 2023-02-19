#------------------------------------------------
# Prototype API for TALOS
#
# Author: Yuan Wang
# Date:  08/08/2022
#------------------------------------------------

home_dir <- "~/SPEEDI"
source(paste0(home_dir, "/prototype_utils.R"))

# Human CC Genes
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Mouse CC Genes
#cc.gene.updated.mouse <- readRDS(paste0(home_dir, "/cc.gene.updated.mouse.rds"))
#m.s.genes <- cc.gene.updated.mouse$m.s.genes
#m.g2m.genes <- cc.gene.updated.mouse$m.g2m.genes

SEED <- 1824409L
set.seed(SEED)

run_SPEEDI <- function(data_path, output_dir, sample_id_list, naming_token, save_progress = TRUE, use_simplified_reference = FALSE, remove_doublets = FALSE) {
  all_sc_exp_matrices <- Read_h5(data_path, sample_id_list)
  assign("sc_obj", FilterRawData(all_sc_exp_matrices, human = TRUE, remove_doublets = FALSE), envir = .GlobalEnv)
  rm(all_sc_exp_matrices)
  assign("sc_obj", InitialProcessing(sc_obj, human = TRUE), envir = .GlobalEnv)
  if(save_progress) {
    save(sc_obj, file = paste0(output_dir, "3_", naming_token, "_sc_obj.rds"))
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
  assign("sc_obj", MapCellTypes(sc_obj, reference), envir = .GlobalEnv)
  rm(reference)
  # Add cell names as a metadata column - this is handy for selecting subsets of cells
  cell_names <- rownames(sc_obj@meta.data)
  assign("sc_obj", AddMetaData(sc_obj, metadata = cell_names, col.name = "cell_name"), envir = .GlobalEnv)
  if(save_progress) {
    saveRDS(sc_obj, file = paste0(output_dir, "7_", naming_token, "_sc_obj.rds"))
  }
  return(sc_obj)
}

print_UMAP <- function(sc_obj, sample_count, group_by_category, output_dir, naming_token, file_suffix) {
  cell_count <- length(sc_obj$cell_name)
  current_title <- paste0("scRNA-seq and/or snRNA-seq Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
  DimPlot(sc_obj, reduction = "umap", group.by = group_by_category, label = TRUE,
        label.size = 3, repel = TRUE, raster = FALSE) +
  labs(title = current_title) +
  theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(output_dir, naming_token, file_suffix), device = "png", dpi = 300)
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
    cluster_prediction <- sort(table(filtered_cluster$predicted_celltype_majority_vote), decreasing = TRUE)[1]
    cluster_predictions <- append(cluster_predictions, cluster_prediction)
    cluster_distributions[[idx]] <- table(filtered_cluster$predicted.id)
    cluster_mean_S_score <- append(cluster_mean_S_score, mean(filtered_cluster$S.Score))
    cluster_mean_G2M_score <- append(cluster_mean_G2M_score, mean(filtered_cluster$G2M.Score))
    cluster_mean_CC_difference <- append(cluster_mean_CC_difference, mean(filtered_cluster$CC.Difference))
    # High
    idxPass <- which(filtered_cluster$viral_load %in% "HIGH")
    cellsPass <- names(filtered_cluster$orig.ident[idxPass])
    num_cells_high <- c(num_cells_high, length(cellsPass))
    if(length(cellsPass) > 0) {
      filtered_cluster_high <- subset(x = filtered_cluster, subset = cell_name %in% cellsPass)
      cluster_mean_mito_high <- c(cluster_mean_mito_high, mean(filtered_cluster_high$percent.mt))
      cluster_mean_nFeature_high <- c(cluster_mean_nFeature_high, mean(filtered_cluster_high$nFeature_RNA))
      cluster_mean_nCount_high <- c(cluster_mean_nCount_high, mean(filtered_cluster_high$nCount_RNA))
      cluster_mean_rp_high <- c(cluster_mean_rp_high, mean(filtered_cluster_high$percent.rp))
    } else {
      cluster_mean_mito_high <- c(cluster_mean_mito_high, NA)
      cluster_mean_nFeature_high <- c(cluster_mean_nFeature_high, NA)
      cluster_mean_nCount_high <- c(cluster_mean_nCount_high, NA)
      cluster_mean_rp_high <- c(cluster_mean_rp_high, NA)
    }
    # Low
    idxPass <- which(filtered_cluster$viral_load %in% "LOW")
    cellsPass <- names(filtered_cluster$orig.ident[idxPass])
    num_cells_low <- c(num_cells_low, length(cellsPass))
    if(length(cellsPass) > 0) {
      filtered_cluster_low <- subset(x = filtered_cluster, subset = cell_name %in% cellsPass)
      cluster_mean_mito_low <- c(cluster_mean_mito_low, mean(filtered_cluster_low$percent.mt))
      cluster_mean_nFeature_low <- c(cluster_mean_nFeature_low, mean(filtered_cluster_low$nFeature_RNA))
      cluster_mean_nCount_low <- c(cluster_mean_nCount_low, mean(filtered_cluster_low$nCount_RNA))
      cluster_mean_rp_low <- c(cluster_mean_rp_low, mean(filtered_cluster_low$percent.rp))
    } else {
      cluster_mean_mito_low <- c(cluster_mean_mito_low, NA)
      cluster_mean_nFeature_low <- c(cluster_mean_nFeature_low, NA)
      cluster_mean_nCount_low <- c(cluster_mean_nCount_low, NA)
      cluster_mean_rp_low <- c(cluster_mean_rp_low, NA)
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

FilterRawData <- function(sc_obj, human, remove_doublets = FALSE) {
  message("Step 2: Filtering out bad samples...")
  testing_flag <- TRUE

  sc_obj <- CreateSeuratObject(counts = all_sc_exp_matrices,
                               assay = "RNA",
                               min.cells = 3,
                               min.features = 3,
                               project = "unbias")

  sc_obj$sample <- as.vector(sapply(strsplit(colnames(sc_obj), "#"), "[", 1))

  if(remove_doublets) {
    message("Removing doublets...")
    sc_obj <- Seurat::as.Seurat(scDblFinder::scDblFinder(Seurat::as.SingleCellExperiment(sc_obj), samples = "sample", BPPARAM=MulticoreParam(7, RNGseed=SEED)))
    # See distribution of doublets in each sample
    doublet_sc_obj <- subset(x = sc_obj, subset = scDblFinder.class %in% "doublet")
    message("Number of doublets removed in each sample:")
    print(table(doublet_sc_obj$sample))
    rm(doublet_sc_obj)
    sc_obj <- subset(x = sc_obj, subset = scDblFinder.class %in% "singlet")
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
    if(!testing_flag) {
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
      object <- subset(objects[[i]], nFeature_RNA > 900 & nFeature_RNA < 4000 & nCount_RNA < 10000 & percent.mt < 15 & percent.hb < 0.4 & percent.rp < 50)
    }
    if(!testing_flag) {
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

#   saveRDS(sc_obj, paste0(home_dir, "/unintegrated.object.RData"))
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
                                      k.anchor = 5)
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




VisualizeIntegration <- function(sc_obj) {
  set.seed(1824409L)
  sc_obj <- ScaleData(sc_obj, verbose = T)
#   sc_obj <- PCA(sc_obj)
  sc_obj <- RunPCA(sc_obj, npcs = 30, approx = T, verbose = T, seed.use = SEED)
  sc_obj <- RunUMAP(sc_obj, reduction = "pca", dims = 1:30, seed.use = SEED, return.model = T)
#   sc_obj <- RunUMAP(sc_obj, reduction = "prcomp", dims = 1:30, seed.use = SEED, return.model = T)
  DefaultAssay(sc_obj) <- "SCT"
  #sc_obj <- PrepSCTFindMarkers(sc_obj)
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
        reference <- LoadH5Seurat(paste0(home_dir, "/reference/pbmc_multimodal.h5seurat"))
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

FindMappingAnchors <- function(sc_obj, reference) {
  DefaultAssay(sc_obj) <- "integrated"
  anchors <- FindTransferAnchors(reference = reference,
                                 query = sc_obj,
                                 normalization.method = "SCT",
                                 recompute.residuals = T,
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
  idx <- grep("cDC", sc_obj$predicted.id)
  sc_obj$predicted.id[idx] <- "cDC"
  #idx <- grep("Proliferating", sc_obj$predicted.id)
  #sc_obj$predicted.id[idx] <- "Proliferating"
  #idx <- grep("B", sc_obj$predicted.id)
  #sc_obj$predicted.id[idx] <- "B"
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "CD4 Naive", "T Naive")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "CD8 Naive", "T Naive")
  sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "NK_CD56bright", "NK")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "ASDC", "CD14 Mono")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "cDC", "CD14 Mono")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "Eryth", "CD14 Mono")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "HSPC", "CD14 Mono")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "pDC", "CD14 Mono")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "Plasmablast", "CD14 Mono")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "Platelet", "CD14 Mono")
  #sc_obj$predicted.id <- replace(sc_obj$predicted.id, sc_obj$predicted.id == "Treg", "T Naive")

  #integrated_snn_res_values <- sc_obj[[associated_res_attribute]]$

  cluster.dump <- as.numeric(levels(sc_obj$integrated_snn_res.1.5))
  sc_obj$predicted_celltype_majority_vote <- sc_obj$seurat_clusters
  levels(sc_obj$predicted_celltype_majority_vote) <- as.character(levels(sc_obj$predicted_celltype_majority_vote))
  for (i in unique(sc_obj$predicted.id)) {
    print(i)
    cells <- names(sc_obj$predicted.id[sc_obj$predicted.id == i])
    freq.table <- as.data.frame(table(sc_obj$integrated_snn_res.1.5[cells]))
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
          cells <- names(sc_obj$integrated_snn_res.1.5[sc_obj$integrated_snn_res.1.5 == i])
          freq.table <- as.data.frame(table(sc_obj$predicted.id[cells]))
          levels(sc_obj$predicted_celltype_majority_vote)[levels(sc_obj$predicted_celltype_majority_vote) %in% as.character(i)] <- as.vector(freq.table$Var1)[which.max(freq.table$Freq)]
      }
  }


  message("...End majority voting")
  return(sc_obj)
}

MapCellTypes <- function(sc_obj, reference) {
  message("Step 6: Reference-based cell type mapping...")
  anchors <- FindMappingAnchors(sc_obj, reference)
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




